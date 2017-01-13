#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>

#include <sys/time.h>


#define true 1
#define false 0
#define HASHSIZE 101

extern int errno;

// Last change: 17-02-2014
// gcc -Wall -Wextra -posix -Werror -std=c99 -D_POSIX_C_SOURCE=2001 -lm  bed_utils_light.c -o bed_utils_light
// valgrind -v --tool=memcheck --leak-check=full ./a.out tiny.bed

struct s_bedinfo {
    int error_code;
    int reads;
    int reads_length_mean;
    int ureads;
    int hreads;
    int sorted;
    char* used_chrm;
    char* last_line;
} bedinfo;

struct s_chrinfo {
    char* chrom;
    unsigned long size;
    unsigned char used;
    unsigned int line_nbr;
    struct s_linedata* head_array;
    struct s_linedata* array_line;
} chrinfo;


struct chrm_entry { 
    char *name; 
    unsigned long size;
    unsigned int used;
    struct chrm_entry *next; 
};


char* strdup(char* s) {
    char* p;
    p = (char*) malloc(strlen(s)+1); 
    if (p != NULL)
       strcpy(p, s);
    return p;
}

// ##################################################################################3

static struct chrm_entry* chrm_hashtab[HASHSIZE]; 
static char **chrm_array;
static int chrm_num;

unsigned int hash(char* s) {
    unsigned int hashval;
    for (hashval = 0; *s != '\0'; s++) {
        hashval = *s + 31 * hashval;
    }
    return hashval % HASHSIZE;
}


struct chrm_entry* lookup_chrm(char* s) {
    struct chrm_entry* ce;
    for (ce = chrm_hashtab[hash(s)]; ce != NULL; ce = ce->next)
        if (strcmp(s, ce->name) == 0)
          return ce; /* found */
    return NULL; /* not found */
}


struct chrm_entry* put_chrm(char* name, unsigned long size) {
    struct chrm_entry *ce;
    unsigned int hashval;
    if ((ce = lookup_chrm(name)) == NULL) { /* not found */
        ce = (struct chrm_entry*)malloc(sizeof(*ce));
        if (ce == NULL) {
            return NULL;
        } else {
            ce->name = strdup(name);
            if (ce->name == NULL) {
                 return NULL;
            }
        }
        hashval = hash(name);
        ce->next = chrm_hashtab[hashval];
        chrm_hashtab[hashval] = ce;
    }
    ce->size = size;
    ce->used = 0;
    if (ce->size == 0) {
       return NULL;
    }
    return ce;
}


void destroy_chrm_entry(struct chrm_entry* cm) {
    free(cm->name);
    free(cm);
}


void remove_chrm(char* name) {
    struct chrm_entry* ce = NULL;
    struct chrm_entry* ceb = NULL;
    for (ce = chrm_hashtab[hash(name)]; ce != NULL; ce = ce->next)
        if (strcmp(name, ce->name) == 0) {
            if (ceb != NULL) {
                ceb->next = ce->next;
            }
            destroy_chrm_entry(ce);
            return;
        }
        ceb = ce;
    return;
}


void destroy_chrm_entry_hash() {
    struct chrm_entry* tmp;
    struct chrm_entry* next;
    for (int i = 0; i < HASHSIZE; i++) {
        tmp = chrm_hashtab[i];
        while(tmp != NULL) {
            next = tmp->next;
            destroy_chrm_entry(tmp);
            tmp = next;
        }
        chrm_hashtab[i] = NULL;
    }
}


//####################################################################################

void get_chrom_size(char* assembly_file) {

    char line[1024];
    char* chrname;
    unsigned char error = 0;
    int i = 0;
    char* token;
    char* ssize;
    unsigned long int size;

    errno = 0;

    FILE* fp = fopen(assembly_file, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: cannot open the assembly file, no such file or directory.\n");
        exit(1);
    }

    char chrm_arr_tmp[100][64];
    
    while(fgets(line, 1024, fp) != NULL)    {
        if (line[0] != '\n' && line[1] != '\n') {
            token = strtok(line,"\t");
            if (token == NULL || strlen(token) == 0) {
                error = 1;
                break;
            }
            chrname = strdup(token);

            token = strtok(NULL,"\t");
            if (token == NULL || strlen(token) == 0) {
                error = 1;
                break;
            }
            size = strtoul(token, &ssize, 10);

            if (errno != 0 || size < 1) {
                error = 1;
                break;
            }
            if (put_chrm(chrname, size) == NULL) {
                error = 1;
                break;
            }
            strcpy(chrm_arr_tmp[i], chrname);
            i++;
        }
    }

    chrm_array = malloc(i * sizeof(char*));
    for (int k = 0; k < i; k++) {
        chrm_array[k] = strdup(chrm_arr_tmp[k]);
    }
    chrm_num = i;

    fclose(fp);

    if (error == 1 || i == 0) {
        fprintf(stderr, "Error: cannot read chromosome information, the assembly file seems to be invalid.\n\
Are you sure the file is tab-delimited?\n");
        exit(1);
    }

    return;
}


struct s_bedinfo* check_bed_file_integrity() {
    char line[1024];
    
    unsigned long start = 0;
    unsigned long end = 0;
    int read_size = 0;
    unsigned long total_read_size = 0;
    unsigned int reads = 0;
    unsigned int header_reads = 0;

    unsigned int reads_length_mean = 0;

    int error_code = 0;
    int line_len = 0;

    char* ptr;
    char* ptrend;

    char strptr_buf[128] = "";

    int lenptr = 0;
    errno = 0;
    char* token;
    int l4p; // true if the is 4p for the current line
    int hl4p = false; //true if one ligne with 4 fields only
    int l5p; // true if the is 5p for the current line
    int hl5p = false; //true if one ligne with 5 fields only

    int sp; //space count
    int vpl; // true if current line has more than 3 space separator
    int hvpl = false; // true if one line une space delimiter

    int vl = false; // true if one ligne is corretly formatted

    struct s_bedinfo* bed_info = (struct s_bedinfo*)malloc(sizeof(struct s_bedinfo));

    int used_chrm_alloc = 256;
    char* used_chrm = (char*)malloc(used_chrm_alloc * sizeof(char));

    struct chrm_entry* ce;

    //error code:
    // 0 the file is okay
    // 1 invalid line(s)
    // 2 strand is in the 4 position
    // 4 strand is in the 5 position
    // 8 space delimiter strand is in the 6 position
    // 16 space delimiter, strand is in the 4 position
    // 32 space delimiter, strand is in the 5 position
    // any other value means different formatted line

    while(fgets(line, 1024, stdin) != NULL) {
        //printf("%s", line);
        line_len = strlen(line);
        //strip the carriage return
        if (line[line_len-1] == '\n') {
            if (line[line_len-2] == '\r') {
                line[line_len-2] = '\0';
            } else {
                line[line_len-1] = '\0';
            }
        }
        if (line[0] == '\n' || line[1] == '\n') {
            continue;
        }
        if (strncmp(line, "#", 1) != 0 && strncmp(line, "track", 5) != 0 && strncmp(line, "browser", 7) != 0) {
            reads++;
            ptr = line;
            vpl = false;
            //parse chrmname -----------------
            if ((ptrend = strchr(ptr, '\t')) == NULL) {
                //no tab character found
                //replace spaces by tab
                sp = 0;
                while ((ptrend = strchr (ptr, ' ')) != NULL) {
                    *ptrend++ = '\t';
                    sp++;
                }
                if (sp < 3) {
                    // no tab delimited and less than 3 space character, abord
                    error_code = 1;
                    break;
                }
                // space have been replace by tab
                vpl = true;
                hvpl = true;
                // get the chrmname again
                ptr = line;
                ptrend = strchr(ptr, '\t');
            }

            //different separator line
            if (hvpl && !vpl) {
                error_code |= 64;
                break;
            }
            lenptr = ptrend-ptr;
            if (lenptr == 0) {
                //empty chrmname, abord
                error_code = 1;
                break;
            }
            strncpy(strptr_buf, ptr, lenptr);
            strptr_buf[lenptr] = '\0';

            //store the chrm used
            if ((ce = lookup_chrm(strptr_buf)) == NULL) {
                put_chrm(strptr_buf, 1); //1 is the size, not used
                while (((int)strlen(strptr_buf) + (int)strlen(used_chrm) + 1) > used_chrm_alloc) {
                    used_chrm_alloc *= 2;
                    if ((used_chrm = realloc(used_chrm, used_chrm_alloc)) == NULL) {
                        fprintf(stderr, "Error: cannot allocate user_chrm memory.\n");
                        exit(1);
                    }
                }
                strcat(strptr_buf, ",");
                strcat(used_chrm, strptr_buf);
            }

            ptr = ptrend+1;
            //parse start position -----------
            if ((ptrend = strchr(ptr, '\t')) == NULL) {
                //no more field, abord
                error_code = 1;
                break;
            }
            lenptr = ptrend-ptr;
            strncpy(strptr_buf, ptr, lenptr);
            strptr_buf[lenptr] = '\0';
            start = strtoul(strptr_buf, &token, 10);
            if (*token != 0 || errno != 0) {
                //invalid numeric value for start position, abord
                error_code = 1;
                break;
            }
            ptr = ptrend+1;
            //parse end position -------------
            if ((ptrend = strchr(ptr, '\t')) == NULL) {
                //no more field, abord
                error_code = 1;
                break;
            }
            lenptr = ptrend-ptr;
            strncpy(strptr_buf, ptr, lenptr);
            strptr_buf[lenptr] = '\0';
            end = strtoul(strptr_buf, &token, 10);
            if (*token != 0 || errno != 0) {
                //invalid numeric value for end position, abord
                error_code = 1;
                break;
            }
            //compute read length
            if ((read_size = (end - start)) < 1) {
                error_code = 1;
                break;
            }
            total_read_size += read_size;
            ptr = ptrend+1;
            // read ID -------------------
            if ((ptrend = strchr(ptr, '\t')) == NULL) {
                //no more field test if strand field
                if (strlen(ptr) == 1 && (*ptr == '+' || *ptr == '-')) {
                    // strand at 4 pos
                    if (vpl) {
                        error_code |= 16;
                    } else {
                        error_code |= 2;
                    }
                    hl4p = true;
                    continue;
                } else {
                    error_code = 1;
                    break;
                }
            }
            //test if the line is longer than other
            l4p = false;
            if ((ptrend-ptr) == 1 && (*ptr == '+' || *ptr == '-')) {
                // strand at pos 4, but keep going parsing until field 6th
                l4p = true;
            }
            if (hl4p && !l4p) {
                error_code |= 64;
                break;
            }
            if (hl4p) {
                continue;
            }
            ptr = ptrend+1;
            //parsing score ------------------
            if ((ptrend = strchr(ptr, '\t')) == NULL) {
                if (strlen(ptr) == 1 && (*ptr == '+' || *ptr == '-')) {
                    // strand at 5 pos
                    if (vpl) {
                        error_code |= 32;
                    } else {
                        error_code |= 4;
                    }
                    hl5p = true;
                } else if (l4p) {
                    if (vpl) {
                        error_code |= 16;
                    } else {
                        error_code |= 2;
                    }
                    hl4p = true;
                } else {
                    error_code = 1;
                    break;
                }
                continue;
            }
            //test if the line is longer than other
            l5p = false;
            if ((ptrend-ptr) == 1 && (*ptr == '+' || *ptr == '-')) {
                // strand at pos 5, but keep going parsing until field 6th
                l5p = true;
            }
            if (hl5p && !l5p) {
                error_code |= 64;
                break;
            }
            if (hl5p) {
                continue;
            }
            ptr = ptrend+1;
            //parse strand -----------------
            if ((ptrend = strchr(ptr, '\t')) == NULL) {
                if (!(strlen(ptr) == 1 && (*ptr == '+' || *ptr == '-'))) {
                    if (l4p) {
                        if (vpl) {
                            error_code |= 32;
                        } else {
                            error_code |= 2;
                        }
                    } else if (l5p) {
                        if (vpl) {
                            error_code |= 16;
                        } else {
                            error_code |= 4;
                        }
                    } else {
                        error_code = 1;
                        break;
                    }
                    continue;
                }
                if (vpl) {
                    if (vl) {
                        error_code |= 64;
                        break;
                    } else {
                        error_code |= 8;
                    }
                } else {
                    vl = true;
                }
            } else {
                if (!((ptrend-ptr) == 1 && (*ptr == '+' || *ptr == '-'))) {
                    error_code = 1;
                    break;
                }
                // valid line
                // get a value for vpl also
                if (vpl) {
                    if (vl) {
                        error_code |= 64;
                        break;
                    } else {
                        error_code |= 8;
                    }
                } else {
                    vl = true;
                }
            }
        } else {
            header_reads++;
        }
    }

    //get reads length mean
    reads_length_mean = total_read_size/reads;

    //fix chrm used string
    if (strlen(used_chrm) > 0) {
        used_chrm[strlen(used_chrm)-1] = '\0';
    } else {
        strcpy(used_chrm, "-");
    }

    bed_info->error_code = error_code;
    bed_info->reads = reads;
    bed_info->reads_length_mean = reads_length_mean;
    bed_info->ureads = 0;
    bed_info->hreads = header_reads;
    bed_info->sorted = false;
    bed_info->used_chrm = used_chrm;
    bed_info->last_line = line;

    return bed_info;
}

//for sort function
int compare(const void *a, const void *b){
    const char **ia = (const char **)a;
    const char **ib = (const char **)b;
    return strcmp(*ia, *ib);
}


//###############################################################3
//dictionnary

struct key_order {
    unsigned int position;
    struct maillon2* ptr;
};

struct maillon2 {
    unsigned int key2;
    unsigned int intensity;
    char* tag;
    struct maillon2* next;
};

struct maillon1* dict1[HASHSIZE];

struct maillon1 {
    char* key1;
    unsigned int nb_keys;
    unsigned int key_array_size;
    struct key_order* key_array;
    struct maillon2** dict2;
    struct maillon1* next;
};

unsigned int hash2(const char* s) {
    unsigned int hashval;
    for (hashval = 0; *s != '\0'; s++) {
        hashval = *s + 31 * hashval;
    }
    return hashval % HASHSIZE;
}

unsigned int hash4(unsigned int key) {
    return key % 1000003;
}


struct maillon1* is_key1(const char* key) {
    struct maillon1* m;
    for (m = dict1[hash2(key)]; m != NULL; m = m->next)
        if (strcmp(m->key1, key) == 0)
            return m;
    return NULL;
}


struct maillon2* is_key2(struct maillon2** dict, const unsigned int key) {
    struct maillon2* m;
    for (m = dict[hash4(key)]; m != NULL; m = m->next) {
        if (m->key2 == key) {
            return m;
        }
    }
    return NULL;
}

struct maillon1* insert_m1(char* key) {
    struct maillon1* m;
    if ((m = is_key1(key)) != NULL) {
        return m;
    } else {
        m = malloc(sizeof(*m));
        if (m == NULL)
            return NULL;
        unsigned int hash = hash2(key);
        m->next = dict1[hash];
        dict1[hash] = m;
        m->key1 = strdup(key);
        m->dict2 = malloc(sizeof(struct maillon2*) * 1000003);
        if (m->dict2 == NULL)
            return NULL;
        m->key_array = malloc(sizeof(struct key_order) * 10000);
        if (m->key_array == NULL)
            return NULL;
        m->nb_keys = 0;
        m->key_array_size = 10000;
    }
    return m;
}

struct maillon2* insert_m2(struct maillon2** dict, const unsigned int key, char* tag, struct maillon1* m1) {
    struct maillon2* m;
    if ((m = is_key2(dict, key)) != NULL) {
        char* new_tag = malloc(strlen(m->tag) + strlen(tag) + 3);
        if (new_tag == NULL) {
            fprintf(stderr, "Error: cannot allocate memory.\n");
            exit(1);
        }
        strcat(new_tag, m->tag);
        strcat(new_tag, ", ");
        strcat(new_tag, tag);
        m->tag = new_tag;
        return m;
    } else {
        m = malloc(sizeof(*m));
        if (m == NULL)
            return NULL;
        unsigned int index = hash4(key);
        m->next = dict[index];
        dict[index] = m;
        m->key2 = key;
        m->intensity = 0;
        m->tag = strdup(tag);
        if (m1->nb_keys == m1->key_array_size) {
            m1->key_array_size = m1->key_array_size*1.25;
            struct key_order* tmp = realloc(m1->key_array, sizeof(struct key_order) * (int)(m1->key_array_size));
            if (tmp == NULL)
                return NULL;
            m1->key_array = tmp;
        }
        (m1->key_array[m1->nb_keys]).position = key;
        (m1->key_array[m1->nb_keys]).ptr = m;
        m1->nb_keys++;
    }
    return m;
}

int sort_array_key(const void *a, const void *b) {
    struct key_order *ia = (struct key_order *)a;
    struct key_order *ib = (struct key_order *)b;
    return ia->position - ib->position;
}

//###############################################################################################################

void bin_coord(FILE* fp_in, unsigned int window_size) {
    char line[1024];
    char chrname[128] = "";
    char transcrit[256] = "";
    char gene[128] = "";

    long start; 
    long end;
    char strand;

    char* ptr;
    char* ptrend;
    char* strtoptr;

    unsigned long start_window;
    unsigned long end_window;

    unsigned int i = 0, count = 0;

    int lenptr;
    char token[128] = "";

    struct maillon1* m1;
    struct maillon2* m2;

    while(fgets(line, 1024, fp_in) != NULL) {
        //printf("%s\n", line);
        if (line[0] == '\n' || line[1] == '\n') {
            continue;
        }
        ptr = line;
        ptrend = strchr(ptr, '\t');
        if (ptrend == NULL) {
            fprintf(stderr, "Error: invalid regions line '%s'.\n", line);
            exit(1);
        }
        lenptr = ptrend-ptr;
        if (lenptr == 0) {
            fprintf(stderr, "Error: invalid regions line '%s'.\n", line);
            exit(1);
        }
        strncpy(chrname, ptr, lenptr);
        chrname[lenptr] = '\0';
        ptr = ptrend+1;

        ptrend = strchr(ptr, '\t');
        if (ptrend == NULL) {
            fprintf(stderr, "Error: invalid regions line '%s'.\n", line);
            exit(1);
        }
        lenptr = ptrend-ptr;
        if (lenptr == 0) {
            fprintf(stderr, "Error: invalid regions line '%s'.\n", line);
            exit(1);
        }
        strncpy(token, ptr, lenptr);
        token[lenptr] = '\0';
        start = strtoul(token, &strtoptr, 10);
        if (*strtoptr != 0 || errno != 0 || start < 0) {
            fprintf(stderr, "Error: invalid regions line '%s'.\n", line);
            exit(1);;
        }
        ptr = ptrend+1;

        ptrend = strchr(ptr, '\t');
        if (ptrend == NULL) {
            fprintf(stderr, "Error: invalid regions line '%s'.\n", line);
            exit(1);
        }
        lenptr = ptrend-ptr;
        if (lenptr == 0) {
            fprintf(stderr, "Error: invalid regions line '%s'.\n", line);
            exit(1);
        }
        strncpy(token, ptr, lenptr);
        token[lenptr] = '\0';
        end = strtoul(token, &strtoptr, 10);
        if (*strtoptr != 0 || errno != 0 || end < 0 || end <= start) {
            fprintf(stderr, "Error: invalid regions line '%s'.\n", line);
            exit(1);;
        }
        ptr = ptrend+1;

        ptrend = strchr(ptr, '\t');
        if (ptrend == NULL) {
            fprintf(stderr, "Error: invalid regions line '%s'.\n", line);
            exit(1);
        }
        lenptr = ptrend-ptr;
        if (lenptr == 0) {
            fprintf(stderr, "Error: invalid regions line '%s'.\n", line);
            exit(1);
        }
        strncpy(transcrit, ptr, lenptr);
        transcrit[lenptr] = '\0';
        ptr = ptrend+1;

        ptrend = strchr(ptr, '\t');
        if (ptrend == NULL) {
            fprintf(stderr, "Error: invalid regions line '%s'.\n", line);
            exit(1);
        }
        lenptr = ptrend-ptr;
        if (lenptr == 0) {
            fprintf(stderr, "Error: invalid regions line '%s'.\n", line);
            exit(1);
        }
        strncpy(gene, ptr, lenptr);
        gene[lenptr] = '\0';
        ptr = ptrend+1;

        ptrend = strchr(ptr, '\t');
        if (ptrend == NULL) {
            if (*ptr == '+') {
                strand = '+';
            } else {
                strand = '-';
            }
            token[0] = strand;
            token[1] = '\0';
        } else {
            lenptr = ptrend-ptr;
            if (lenptr == 0) {
                fprintf(stderr, "Error: invalid regions line '%s'.\n", line);
                exit(1);
            }
            strncpy(token, ptr, lenptr);
            token[lenptr] = '\0';
            if (token[0] == '+') {
                strand = '+';
            } else {
                strand = '-';
            }
        }

        //bin the regions
        strcat(transcrit,":");
        strcat(transcrit, gene);
        strcat(transcrit,":");
        strcat(transcrit, token);

        if ((m1 = is_key1(chrname)) == NULL) {
            m1 = insert_m1(chrname);
            if (m1 == NULL) {
                fprintf(stderr, "Error: cannot allocate memory block.\n");
                exit(1);
            }
        }


        start_window = (start / window_size) * window_size;
        end_window = (end / window_size) * window_size;
        for (i = start_window; i <= end_window; i+= window_size) {
            m2 = insert_m2(m1->dict2, i, transcrit, m1);
            count++;
            if (m2 == NULL) {
                fprintf(stderr, "Error: cannot allocate memory block.\n");
                exit(1);
            }
        }
    }
}


void pileup(FILE* fp_in, FILE* fp_out, unsigned int window_size, unsigned int extend_size, unsigned char no_header, unsigned char remove_clonal, unsigned char middle_pos) {

    char line[1024];
    char chrname[128] = "";
    char current_chrm[128] = "";
    unsigned long start_pos; 
    long tmp;
    unsigned long end_pos;
    char strand;

    char* ptr;
    char* ptrend;

    unsigned long start_window;
    unsigned long end_window;


    unsigned int prev_start_pos = 0;
    unsigned int prev_end_pos = 0;
    char prev_strand = 0;


    unsigned int read_length = 0;
    int distance = 0;

    unsigned int r = 0, i = 0;

    int lenptr;
    char token[128] = "";

    struct maillon1* m1, *m1_tmp;
    struct maillon2* m2;

    while(fgets(line, 1024, stdin) != NULL) {
        if (line[0] == '\n' || line[1] == '\n') {
            continue;
        } else if (no_header || (strncmp(line, "#", 1) != 0 && strncmp(line, "track", 5) != 0 && strncmp(line, "browser", 7) != 0)) {
            //printf("%s\n", line);
            ptr = line;
            ptrend = strchr(ptr, '\t');
            lenptr = ptrend-ptr;
            strncpy(chrname, ptr, lenptr);
            chrname[lenptr] = '\0';
            ptr = ptrend+1;

            if ((m1 = is_key1(chrname)) == NULL) {
                continue;
            }

            ptrend = strchr(ptr, '\t');
            lenptr = ptrend-ptr;
            strncpy(token, ptr, lenptr);
            token[lenptr] = '\0';
            start_pos = atol(token);
            ptr = ptrend+1;

            ptrend = strchr(ptr, '\t');
            lenptr = ptrend-ptr;
            strncpy(token, ptr, lenptr);
            token[lenptr] = '\0';
            end_pos = atol(token) - 1;

            if (extend_size != 0) {
                read_length = end_pos - start_pos;
                distance = extend_size - read_length;

                if (distance > 0) {
                    ptr = ptrend+1;
                    ptrend = strchr(ptr, '\t');
                    ptr = ptrend+1;
                    ptrend = strchr(ptr, '\t');
                    ptr = ptrend+1;
                    ptrend = strchr(ptr, '\t');
                    if (ptrend == NULL) {
                        if (*ptr == '+') {
                            strand = '+';
                        } else {
                            strand = '-';
                        }
                    } else {
                        lenptr = ptrend-ptr;
                        strncpy(token, ptr, lenptr);
                        token[lenptr] = '\0';
                        if (token[0] == '+') {
                            strand = '+';
                        } else {
                            strand = '-';
                        }
                    }

                    if (strand == '+') {
                        end_pos = end_pos + distance;
                    } else {
                        tmp = start_pos - distance;
                        if (tmp < 0) {
                            start_pos = 0; 
                        } else {
                            start_pos = tmp;
                        }
                    }

                } else if (remove_clonal) {
                    ptr = ptrend+1;
                    ptrend = strchr(ptr, '\t');
                    ptr = ptrend+1;
                    ptrend = strchr(ptr, '\t');
                    ptr = ptrend+1;
                    ptrend = strchr(ptr, '\t');
                    if (ptrend == NULL) {
                        if (*ptr == '+') {
                            strand = '+';
                        } else {
                            strand = '-';
                        }
                    } else {
                        lenptr = ptrend-ptr;
                        strncpy(token, ptr, lenptr);
                        token[lenptr] = '\0';
                        if (token[0] == '+') {
                            strand = '+';
                        } else {
                            strand = '-';
                        }
                    }
                }
            } else if (remove_clonal) {
                ptr = ptrend+1;
                ptrend = strchr(ptr, '\t');
                ptr = ptrend+1;
                ptrend = strchr(ptr, '\t');
                ptr = ptrend+1;
                ptrend = strchr(ptr, '\t');
                if (ptrend == NULL) {
                    if (*ptr == '+') {
                        strand = '+';
                    } else {
                        strand = '-';
                    }
                } else {
                    lenptr = ptrend-ptr;
                    strncpy(token, ptr, lenptr);
                    token[lenptr] = '\0';
                    if (token[0] == '+') {
                        strand = '+';
                    } else {
                        strand = '-';
                    }
                }
            }

            if (strcmp(chrname, current_chrm) == 0) {

                if (remove_clonal && prev_start_pos == start_pos && prev_end_pos == end_pos && prev_strand == strand) {
                    continue;
                }

                if (middle_pos) {
                    start_window = ((start_pos + (end_pos - start_pos) / 2) / window_size) * window_size;
                    end_window = start_window;
                } else {
                    start_window = (start_pos / window_size) * window_size;
                    end_window = (end_pos / window_size) * window_size;
                }
                for (r = start_window; r <= end_window; r++) {
                    if ((m2 = is_key2(m1->dict2, r)) != NULL) {
                        m2->intensity++;
                    }
                }
            } else {
                if (strncmp(current_chrm, "", 1) != 0) {
                    qsort(m1_tmp->key_array, sizeof(m1_tmp->key_array) / sizeof(struct key_order), sizeof(struct key_order), sort_array_key);
                    //printf("chrm print : %s\n", current_chrm);
                    for (i = 0; i < m1_tmp->nb_keys; i++) {
                        //printf("i : %d\n", i);
                        m2 = (m1_tmp->key_array[i]).ptr;
                        if (m2 != NULL) {
                            fprintf(fp_out, "%s\t%d\t%d\t%d\t%s\n", current_chrm, m2->key2, m2->key2 + window_size, m2->intensity, m2->tag);
                        }
                    }
                }

                strcpy(current_chrm, chrname);

                //the first read
                if (middle_pos) {
                    start_window = ((start_pos + (end_pos - start_pos) / 2) / window_size) * window_size;
                    end_window = start_window;
                } else {
                    start_window = (start_pos / window_size) * window_size;
                    end_window = (end_pos / window_size) * window_size;
                }
                for (r = start_window; r <= end_window; r++) {
                    if ((m2 = is_key2(m1->dict2, r)) != NULL) {
                        m2->intensity++;
                    }
                }
            }
            m1_tmp = m1;
            if (remove_clonal) {
                prev_start_pos = start_pos;
                prev_end_pos = end_pos;
                prev_strand = strand;
            }
        }
    }

    if (strncmp(current_chrm, "", 1) != 0) {
        qsort(m1->key_array, sizeof(m1->key_array) / sizeof(struct key_order), sizeof(struct key_order), sort_array_key);
        for (i = 0; i < m1_tmp->nb_keys; i++) {
            //printf("i : %d\n", i);
            m2 = (m1_tmp->key_array[i]).ptr;
            if (m2 != NULL) {
                fprintf(fp_out, "%s\t%d\t%d\t%d\t%s\n", current_chrm, m2->key2, m2->key2 + window_size, m2->intensity, m2->tag);
            }
        }
    }
fclose(fp_in);
/*fclose(fp_out);*/
    return;
}

void bed_to_wig(FILE* fp, FILE* fp_out, int window_size, int elongation, unsigned char no_header, unsigned char remove_clonal, unsigned char middle_pos) {

    char line[1024];
    char chrname[128] = "";
    char current_chrm[128] = "";
    unsigned long start_pos;
    long tmp;
    unsigned long end_pos;
    char strand;

    char* ptr;
    char* ptrend;

    unsigned long last_index;
    unsigned long last_index_fb;
    unsigned long start_window;
    unsigned long end_window;
    unsigned long start_pos_window;
    unsigned long current_chrm_size = 0;
    unsigned long current_chrm_size_fb = 0;

    unsigned int prev_start_pos = 0;
    unsigned int prev_end_pos = 0;
    char prev_strand = 0;

    unsigned int read_length = 0;
    int distance = 0;

    unsigned int r = 0;

    int lenptr;
    char token[128] = "";

    int* current_array_int = NULL;
    int chrm_array_index = 0;

    qsort(chrm_array, chrm_num, sizeof(char *), compare);

    while(fgets(line, 1024, stdin) != NULL) {
        if (line[0] == '\n' || line[1] == '\n') {
            continue;
        }
        if (no_header || (strncmp(line, "#", 1) != 0 && strncmp(line, "track", 5) != 0 && strncmp(line, "browser", 7) != 0)) {
            //printf("%s\n", line);
            ptr = line;
            ptrend = strchr(ptr, '\t');
            lenptr = ptrend-ptr;
            strncpy(chrname, ptr, lenptr);
            chrname[lenptr] = '\0';
            ptr = ptrend+1;

            if (lookup_chrm(chrname) == NULL) {
                continue;
            }

            ptrend = strchr(ptr, '\t');
            lenptr = ptrend-ptr;
            strncpy(token, ptr, lenptr);
            token[lenptr] = '\0';
            start_pos = atol(token);
            ptr = ptrend+1;

            ptrend = strchr(ptr, '\t');
            lenptr = ptrend-ptr;
            strncpy(token, ptr, lenptr);
            token[lenptr] = '\0';
            end_pos = atol(token) - 1;

            if (elongation != 0) {
                read_length = end_pos - start_pos;
                distance = elongation - read_length;
                if (distance > 0) {
                    ptr = ptrend+1;
                    ptrend = strchr(ptr, '\t');
                    ptr = ptrend+1;
                    ptrend = strchr(ptr, '\t');
                    ptr = ptrend+1;
                    ptrend = strchr(ptr, '\t');
                    if (ptrend == NULL) {
                        if (*ptr == '+') {
                            strand = '+';
                        } else {
                            strand = '-';
                        }
                    } else {
                        lenptr = ptrend-ptr;
                        strncpy(token, ptr, lenptr);
                        token[lenptr] = '\0';
                        if (token[0] == '+') {
                            strand = '+';
                        } else {
                            strand = '-';
                        }
                    }

                    if (strand == '+') {
                        end_pos = end_pos + distance;
                    } else {
                        tmp = start_pos - distance;
                        if (tmp < 0) {
                            start_pos = 0; 
                        } else {
                            start_pos = tmp;
                        }
                    }

                } else if (remove_clonal) {
                    ptr = ptrend+1;
                    ptrend = strchr(ptr, '\t');
                    ptr = ptrend+1;
                    ptrend = strchr(ptr, '\t');
                    ptr = ptrend+1;
                    ptrend = strchr(ptr, '\t');
                    if (ptrend == NULL) {
                        if (*ptr == '+') {
                            strand = '+';
                        } else {
                            strand = '-';
                        }
                    } else {
                        lenptr = ptrend-ptr;
                        strncpy(token, ptr, lenptr);
                        token[lenptr] = '\0';
                        if (token[0] == '+') {
                            strand = '+';
                        } else {
                            strand = '-';
                        }
                    }
                }
            } else if (remove_clonal) {
                ptr = ptrend+1;
                ptrend = strchr(ptr, '\t');
                ptr = ptrend+1;
                ptrend = strchr(ptr, '\t');
                ptr = ptrend+1;
                ptrend = strchr(ptr, '\t');
                if (ptrend == NULL) {
                    if (*ptr == '+') {
                        strand = '+';
                    } else {
                        strand = '-';
                    }
                } else {
                    lenptr = ptrend-ptr;
                    strncpy(token, ptr, lenptr);
                    token[lenptr] = '\0';
                    if (token[0] == '+') {
                        strand = '+';
                    } else {
                        strand = '-';
                    }
                }
            }

            if (strcmp(chrname, current_chrm) == 0) {

                if (remove_clonal && prev_start_pos == start_pos && prev_end_pos == end_pos && prev_strand == strand) {
                    continue;
                }

                if (middle_pos) {
                    start_window = (start_pos + (end_pos - start_pos) / 2) / window_size;
                    end_window = start_window;
                } else {
                    start_window = start_pos / window_size;
                    end_window = end_pos / window_size;
                    if (end_window > last_index) {
                        end_window = last_index;
                    }
                }

                if (start_window > last_index) {
                    fprintf(stderr, "Warning: Coordinate out of bound: %s\t%ld.\nSkipped..", chrname, start_pos);
                    start_window = end_window = last_index;
/*                    if (current_array_int != NULL) {*/
/*                        free(current_array_int);*/
/*                    }*/
/*                    fclose(fp);*/
/*                    fclose(fp_out);*/
/*                    exit(1);*/
                }
                for (r = start_window; r <= end_window; r++) {
                    current_array_int[r] += 1;
                }
            } else {

                while (chrm_array_index < chrm_num && strcmp(chrname, chrm_array[chrm_array_index]) != 0) {
                    if (strcmp(current_chrm, chrm_array[chrm_array_index]) == 0) {
                        ++chrm_array_index;
                        break;
                    }
                    current_chrm_size_fb = lookup_chrm(chrm_array[chrm_array_index])->size;
                    last_index_fb = current_chrm_size_fb / window_size;
                    for (unsigned int i = 0; i < last_index_fb + 1; i++) {
                        start_pos_window = i * window_size;
/*                        end_pos_window = (i * window_size) + window_size;*/
                        fprintf(fp_out, "%s\t%ld\t%ld\t0\n", chrm_array[chrm_array_index], start_pos_window, start_pos_window + window_size-1);
                    }
                    ++chrm_array_index;
                }

                if (strncmp(current_chrm, "", 1) != 0) {
                    for (unsigned int i = 0; i < last_index + 1; i++) {
                        int v = current_array_int[i];
                        start_pos_window = i * window_size;

                        fprintf(fp_out, "%s\t%ld\t%ld\t%d\n", current_chrm, start_pos_window, start_pos_window + window_size-1, v);
                    }
                }

                strcpy(current_chrm, chrname);
                current_chrm_size = lookup_chrm(chrname)->size;
                last_index = current_chrm_size / window_size;
                if (current_array_int != NULL) {
                    free(current_array_int);
                }

                current_array_int = calloc((last_index + 1), sizeof(int));
                if (current_array_int == NULL) {
                    fprintf(stderr, "Error: cannot allocate memory block for current_array_int.\n");
                    fclose(fp);
                    fclose(fp_out);
                    exit(1);
                }

                if (middle_pos) {
                    start_window = (start_pos + (end_pos - start_pos) / 2) / window_size;
                    end_window = start_window;
                } else {
                    start_window = start_pos / window_size;
                    end_window = end_pos / window_size;
                    if (end_window > last_index) {
                        end_window = last_index;
                    }
                }

                if (start_window > last_index) {
                    fprintf(stderr, "Warning: Coordinate out of bound: %s\t%ld.\nSkipped..", chrname, start_pos);
                    start_window = end_window = last_index;
/*                    if (current_array_int != NULL) {*/
/*                        free(current_array_int);*/
/*                    }*/
/*                    fclose(fp);*/
/*                    fclose(fp_out);*/
/*                    exit(1);*/
                }
                for (r = start_window; r <= end_window; r++) {
                    current_array_int[r] += 1;
                }
            }
            if (remove_clonal) {
                prev_start_pos = start_pos;
/*                prev_k*/
                end_pos = end_pos;
                prev_strand = strand;
            }
        }
    }

    //write the last chrmosome info if any


    while (chrm_array_index < chrm_num && strcmp(chrname, chrm_array[chrm_array_index]) != 0) {
        if (strcmp(current_chrm, chrm_array[chrm_array_index]) == 0) {
            ++chrm_array_index;
            break;
        }
        current_chrm_size_fb = lookup_chrm(chrm_array[chrm_array_index])->size;
        last_index_fb = current_chrm_size_fb / window_size;
        for (unsigned int i = 0; i < last_index_fb + 1; i++) {
            start_pos_window = i * window_size;
            fprintf(fp_out, "%s\t%ld\t%ld\t0\n", chrm_array[chrm_array_index], start_pos_window, start_pos_window + window_size-1);
        }
        ++chrm_array_index;
    }

    if (strncmp(current_chrm, "", 1) != 0) {
        for (unsigned int i = 0; i < last_index + 1; i++) {
            int v = current_array_int[i];
            start_pos_window = i * window_size;

            fprintf(fp_out, "%s\t%ld\t%ld\t%d\n", current_chrm, start_pos_window, start_pos_window + window_size-1, v);
        }
        ++chrm_array_index;
    }


    while (chrm_array_index < chrm_num) {
        current_chrm_size_fb = lookup_chrm(chrm_array[chrm_array_index])->size;
        last_index_fb = current_chrm_size_fb / window_size;
        for (unsigned int i = 0; i < last_index_fb + 1; i++) {
            start_pos_window = i * window_size;
            fprintf(fp_out, "%s\t%ld\t%ld\t0\n", chrm_array[chrm_array_index], start_pos_window, start_pos_window + window_size-1);
        }
        ++chrm_array_index;
    }

    fclose(fp);
    fclose(fp_out);
    if (current_array_int != NULL) {
        free(current_array_int);
    }

    return;
}

void check_bed_file_sorted(unsigned char no_header) {
    char line[1024];
    char chrname[32] = "";
    char tmp_chrname[32] = "";
    unsigned long start;
    unsigned long tmp_start = -999999999;
    unsigned char strand;
    unsigned char tmp_strand;

    unsigned int reads = 0;

    char* ptr;
    char* ptrend;

    int cs = 0;
    int lenptr;

    char token[32] = "";

    unsigned int dif_strand = false;
    unsigned int sorted = 1;

    struct chrm_entry* ce;

    while(fgets(line, 1024, stdin) != NULL) {
        //printf("%s",line);
        if (line[0] == '\n' || line[1] == '\n') {
            continue;
        } else if (no_header || (strncmp(line, "#", 1) != 0 && strncmp(line, "track", 5) != 0 && strncmp(line, "browser", 7) != 0)) {
            ptr = line;
            ptrend = strchr(ptr, '\t');
            lenptr = ptrend-ptr;
            strncpy(chrname, ptr, lenptr);
            chrname[lenptr] = '\0';
            ptr = ptrend+1;

            ptrend = strchr(ptr, '\t');
            lenptr = ptrend-ptr;
            strncpy(token, ptr, lenptr);
            token[lenptr] = '\0';
            start = atol(token);

            cs = strcmp(chrname, tmp_chrname);
            //if (cs < 0) { /* chrm not sorted (may not work it sorted base on chrm number, chr9 > chr10 > chr1) */
            //    sorted = 0;
            //    break; 
            /*} else*/ if (cs != 0) { /* different chr */
                if ((ce = lookup_chrm(chrname)) != NULL) { /* if the chrm already exists in the dict, file is unsorted */
                    sorted = 0;
                    break;
                } else {
                     ce = put_chrm(chrname, 1);
                }
            } else {
                if (start < tmp_start) {
                    sorted = 0;
                    break;
                } else if (reads == 0) {
                    ce = put_chrm(chrname, 1); /* first line case*/
                } else if (start == tmp_start) {
                    //check if the strand is different
                    //get the strand
                    ptr = ptrend+1;
                    ptrend = strchr(ptr, '\t');
                    ptr = ptrend+1;
                    ptrend = strchr(ptr, '\t');
                    ptr = ptrend+1;
                    ptrend = strchr(ptr, '\t');
                    ptr = ptrend+1;
                    ptrend = strchr(ptr, '\t');
                    if (ptrend == NULL) {
                        if (*ptr == '+') {
                            strand = 0;
                        } else {
                            strand = 1;
                        }
                    } else {
                        lenptr = ptrend-ptr;
                        strncpy(token, ptr, lenptr);
                        //token[lenptr] = '\0';
                        if (token[0] == '+') {
                            strand = 0;
                        } else {
                            strand = 1;
                        }
                    }
                    if (tmp_strand != strand) {
                        if (dif_strand) {
                            /* if the strand is different and if was different before e.g:
                             * chr1 100 200 tag 0 +
                             * chr1 100 200 tag 0 -
                             * chr1 100 200 tag 0 +
                             * the file is not sorted
                             */
                            sorted = 0;
                            break;
                        }
                        dif_strand = true;
                    }
                    strcpy(tmp_chrname, chrname);
                    tmp_start = start;
                    tmp_strand = strand;
                    reads++;
                    continue;
                } else {
                    /* reset */
                    dif_strand = false;
                }

            }
            //get the strand
            ptr = ptrend+1;
            ptrend = strchr(ptr, '\t');
            ptr = ptrend+1;
            ptrend = strchr(ptr, '\t');
            ptr = ptrend+1;
            ptrend = strchr(ptr, '\t');
            ptr = ptrend+1;
            ptrend = strchr(ptr, '\t');
            if (ptrend == NULL) {
                if (*ptr == '+') {
                    strand = 0;
                } else {
                    strand = 1;
                }
            } else {
                lenptr = ptrend-ptr;
                strncpy(token, ptr, lenptr);
                //token[lenptr] = '\0';
                if (token[0] == '+') {
                    strand = 0;
                } else {
                    strand = 1;
                }
            }
            strcpy(tmp_chrname, chrname);
            tmp_start = start;
            tmp_strand = strand;
            reads++;
        }
    }

    printf("%d\n", sorted);

    return;
}

int main(int argc, char* argv[]) {
    unsigned int uargc = argc;
    //int retCode;
    unsigned int i = 0;
    char *endptr;

    while (i < uargc) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            fprintf(stderr, "usage: %s [integrity|wig|countbin] bed_file [arguments] [options]\n\n\
\tintegrity:\n\nusage: %s integrity bed_file\noutput: [error code] [read count] [mean reads length]\
 [header line count] [last line]\n\n\twig:\n\nusage: %s\
 wig bed_file output_file genome_file window_size elongation [--noheader]\noutput: none\n\n", argv[0], argv[0], argv[0]);
            exit(0);
        }
        i++;
    }

    FILE* fp;

    if (uargc < 3) {
        fprintf(stderr, "usage: %s [integrity|wig|countbin] bed_file [arguments] [options]\n\n\
\tintegrity:\n\nusage: %s integrity bed_file\n\n\twig:\n\nusage: %s wig bed_file output_file genome_file window_size elongation\n\n", argv[0], argv[0], argv[0]);
        exit(1);
    } else {

        fp = fopen(argv[2], "r");
        if (fp == NULL) {
            fputs("Error: cannot open the input file, no such file or directory.\n", stderr);
            exit(1);
        }

        if (strcmp(argv[1], "integrity") == 0) {
            struct s_bedinfo* bi = check_bed_file_integrity();
            printf("%d,%d %d %d,%s\n", bi->error_code, bi->reads, bi->reads_length_mean, bi->hreads, bi->last_line);
            free(bi->used_chrm);
            fclose(fp);
        } else if (strcmp(argv[1], "wig") == 0) {
            unsigned char no_header = false;
            if (argc >= 9 && strcmp(argv[8], "--noheader") == 0) {
                no_header = true;
            }

            unsigned char remove_clonal = false;
            if (argc >= 10 && strcmp(argv[9], "--noclonal") == 0) {
                remove_clonal = true;
            }

            unsigned char middle_pos = false;
            if (argc == 11 && strcmp(argv[10], "--middlepos") == 0) {
                middle_pos = true;
            }


            get_chrom_size(argv[4]);

            // add error check
            int window_size = strtol(argv[5], &endptr, 10);
            int elongation = strtol(argv[6], &endptr, 10);

            FILE* fp_out = fopen(argv[3], "w");
            if (fp_out == NULL) {
                printf("Error: cannot create the output file.\n");
            }
            fprintf(fp_out, "track type=bedGraph name=%s\n", argv[7]);
            bed_to_wig(fp, fp_out, window_size, elongation, no_header, remove_clonal, middle_pos);
        } else if (strcmp(argv[1], "issorted") == 0) {
            unsigned int no_header = false;
            if (argc >= 9 && strcmp(argv[8], "--noheader") == 0) {
                no_header = true;
            }
            check_bed_file_sorted(no_header);
        } else if (strcmp(argv[1], "countbin") == 0) {
            FILE* fp2 = fopen(argv[4], "r");
            if (fp2 == NULL) {
                fputs("Error: cannot open the coordinate file, no such file or directory.\n", stderr);
                exit(1);
            }

            FILE* fp_out = fopen(argv[3], "w");
            if (fp_out == NULL) {
                printf("Error: cannot create the output file.\n");
            }

            unsigned int window_size = strtol(argv[5], &endptr, 10);
            unsigned int elongation = strtol(argv[6], &endptr, 10);

            bin_coord(fp2, window_size);
            fclose(fp2);
            unsigned char no_header = false;
            if (argc >= 8 && strcmp(argv[7], "--noheader") == 0) {
                no_header = true;
            }

            unsigned char remove_clonal = false;
            if (argc >= 9 && strcmp(argv[8], "--noclonal") == 0) {
                remove_clonal = true;
            }

            unsigned char middle_pos = false;
            if (argc == 10 && strcmp(argv[9], "--middlepos") == 0) {
                middle_pos = true;
            }

            pileup(fp, fp_out, window_size, elongation, no_header, remove_clonal, middle_pos);

        } else {
            fputs("Error: invalid task.\n", stderr);
            exit(1);
        }
    }
    return (EXIT_SUCCESS);
}

