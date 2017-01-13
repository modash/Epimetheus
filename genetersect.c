/* 
 * File:   genetersect.c
 * Author: blum
 *
 * Created on 2 février 2015, 19:22
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "uthash.h"


typedef struct struct_gene {
    unsigned long pos;
    char name[1024];
} struct_gene;


typedef struct array {
    struct_gene *elements;
    size_t used;
    size_t size;
} array;

typedef struct dict {
    char chrm[512];
    array a;
    UT_hash_handle hh;
} dict;

void rtrim(char *str) {
    size_t n;
    n = strlen(str);
    while (n > 0 && isspace((unsigned char) str[n - 1])) {
        n--;
    }
    str[n] = '\0';
}


void init(array *a, size_t maxsize) {
    a->elements = malloc(maxsize * sizeof(struct_gene));
    
    if (a->elements == NULL) {
        exit(1);
    }
    a->used = 0;
    a->size = maxsize;
}

void append(array *a, struct_gene *gene) {
    if (a->used == a->size) {
        a->size *= 2;
        a->elements = (struct_gene *)realloc(a->elements, a->size * sizeof(struct_gene));
        if (a->elements == NULL) {
            exit(1);
        }
    }
    a->elements[a->used++] = *gene;
}

void empty(array *a) {
    free(a->elements);
    a->elements = NULL;
    a->used = 0;
    a->size = 0;
}

void intersect(FILE* fp, dict* d) {
    char buf[1024];
    char *str, *token, *ptr1;
    int i;
    unsigned long c;
    int lock = 0;
    
    char *chrm;
    unsigned long pos1, pos2;
    
    while (fgets(buf, sizeof(buf), fp)) {
        i = 0;
        if (strncmp(buf, "track", 5) == 0 && strncmp(buf, "bedGraph", 7) == 0) {
            continue;
            }
        for (str = buf; ; i++) {
            token = strsep(&str, "\t");
            if (token == NULL) {
                break;
            } else if (i == 0) {
                chrm = token;
            } else if (i == 1) {
                pos1 = strtol(token, &ptr1, 10);
                if (ptr1 == token) {
                    exit(1);
                }
            } else if (i == 2) {
                pos2 = strtol(token, &ptr1, 10);
                if (ptr1 == token) {
                    exit(1);
                }                
            } else if (i == 3) {
                if (strcmp(d->chrm, chrm) == 0) {
                    if (! lock) {
                        lock++;
                    }
                    
                    rtrim(token);
                    //printf("-> %s\t%ld\t%ld\t%s\n", chrm, pos1, pos2, token);
                    
                    for (c = 0; c < d->a.used; c++) {
                        if (pos1 <= d->a.elements[c].pos && pos2 >= d->a.elements[c].pos) {
                            printf("%s\t%ld\t%ld\t%s\n", chrm, pos1, pos2, d->a.elements[c].name);
                        }
                    }
                } else if (lock) {
                }

                break;
            }
        }
    }    
}

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "\nUsage: intersect <input1> <input2>\n");
        exit(1);
    }
    
    FILE* fp1 = fopen(argv[1], "r");
    if (fp1 == NULL) {
        fprintf(stderr, "Error: invalid input file (%s)\n", argv[1]);
        exit(1);
    }
    
    FILE* fp2 = fopen(argv[2], "r");
    if (fp2 == NULL) {
        fprintf(stderr, "Error: invalid input file (%s)\n", argv[2]);
        exit(1);
    }
    
    char buf[1024];
    char *str, *str2, *token, *ptr1;
    int integrity = 1;
    int i;
    
    char *chrm, *prev_chrm = NULL, *gene_name;
    unsigned long pos1, pos2;
    
    dict *chroms = NULL;
    dict *d, *o;
    
    while (fgets(buf, sizeof(buf), stdin)) {
        i = 0;
        str2 = strdup(buf);
        rtrim(str2);
        
        for (str = buf; ; i++) {
            token = strsep(&str, "\t");
            
            if (token == NULL) {
                fprintf(stderr, "Illformated line: %s\n", str2);
                integrity = 0;
                break;
            } else if (i == 0) {
                chrm = token;
            } else if (i == 1) {
                pos1 = strtol(token, &ptr1, 10);
                if (ptr1 == token) {
                    exit(1);
                }
            } else if (i == 2) {
                pos2 = strtol(token, &ptr1, 10);
                if (ptr1 == token) {
                    exit(1);
                }                
            } else if (i == 3) {
                gene_name = token;
                rtrim(gene_name);
                
                struct_gene *gene = malloc(sizeof(struct_gene));
                strcpy(gene->name, gene_name);
                gene->pos = pos1;
                
                HASH_FIND_STR(chroms, chrm, d);
                if (! d) {
                    if (prev_chrm != NULL) {
                        HASH_FIND_STR(chroms, prev_chrm, o);
                        intersect(fp2, o);
                        rewind(fp2);
                        HASH_DEL(chroms, o);
                        free(o);
                    }

                    dict *newchrom = malloc(sizeof(dict));
                    strcpy(newchrom->chrm, chrm);
                    init(&newchrom->a, 1000);
                    append(&newchrom->a, gene);
                    HASH_ADD_STR(chroms, chrm, newchrom);
                } else {
                    append(&d->a, gene);
                }
                
                prev_chrm = strdup(chrm);
                break;
            }
        }
        
        if (! integrity) {
            break;
        }
    }
    
    HASH_FIND_STR(chroms, prev_chrm, o);
    if (o) {
        intersect(fp2, o);
        HASH_DEL(chroms, o);
        free(o);            
    }
    
    fclose(fp1);
    fclose(fp2);
    
    return (EXIT_SUCCESS);
}

