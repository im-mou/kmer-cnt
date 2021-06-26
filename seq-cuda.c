
#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fcntl.h>
#include <limits.h>
#include <string.h>

typedef struct __ReadSeqList {
	char* sequence;
	unsigned length;
	struct __ReadSeqList* next;
} ReadSeqList;

typedef struct HashTable {
	unsigned int bits;
	unsigned int count;
	unsigned int read_count;
	//unsigned int *collition;
	unsigned long long int *keys;
    unsigned int *values;
} HashTable;



const unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

HashTable* HashTable_init(unsigned int bits, unsigned int read_count){
	unsigned int capacity = 1U << bits;
    HashTable *ht;
    ht = (HashTable*)calloc(1, sizeof(HashTable));
    ht->keys = (unsigned long long int*)calloc(capacity, sizeof(unsigned long long int));
    ht->values = (unsigned int*)calloc(capacity, sizeof(unsigned int));

	ht->bits = bits;
	ht->count = 0;

    return ht;
}


void HashTable_destory(HashTable *ht) {
	if (!ht) return;
	free((void *)ht->keys);
	free((void *)ht->values);
	free(ht);
}


// funcion para calcular un hash de 64 bits
unsigned int hash_uint64(unsigned long long int key) {

	key = ~key + (key << 21);
	key = key ^ key >> 24;
	key = (key + (key << 3)) + (key << 8);
	key = key ^ key >> 14;
	key = (key + (key << 2)) + (key << 4);
	key = key ^ key >> 28;
	key = key + (key << 31);
	return (unsigned int)key;
}

static inline unsigned int h2b(unsigned int hash, unsigned int bits) {
    return hash * 2654435769U >> (32 - bits);
}

void hash_insert(HashTable *ht, unsigned long long int kmer) {

	unsigned int iKey, last;
	//bool end = false;

	//iKey = last = hash_uint64(kmer) * ((2654435769U >> (32 - ht->bits)));
    iKey = last = h2b(hash_uint64(kmer), ht->bits);

    while (ht->values[iKey] > 0 && ht->keys[iKey] != kmer) {
        iKey = (iKey + 1U) & ((1U << ht->bits) - 1);
        if (iKey == last) break;
    }
    // Comprobar si se ha encontrado un slot vacío
    if (ht->values[iKey] == 0) { // no se ha encontrado la llave

        ht->keys[iKey] = kmer;
        ht->values[iKey] = 1;
        ++ht->count;

    } else {
        ht->values[iKey]++;
    } 
}

// insert k-mers in $seq to hash table $ht
void kernel_count_seq_kmers(HashTable *ht, int k, char *seq)
{

    int i, l;
    //int len = strlen(seq);
    int len = 100;
    unsigned long long int x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;

    for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
        int c = seq_nt4_table[(unsigned char)seq[i]];
        if (c < 4) { // not an "N" base
            x[0] = (x[0] << 2 | c) & mask;                  // forward strand
            x[1] = x[1] >> 2 | (unsigned long long int)(3 - c) << shift;  // reverse strand
            if (++l >= k) { // we find a k-mer

                unsigned long long int kmer = x[0] < x[1]? x[0] : x[1];
                hash_insert(ht, kmer); // only add one strand!

            }
        } else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
    }
}

void kernel_print_hist(const HashTable *ht)
{
    uint32_t k;
    uint64_t cnt[256];
    int i;
    for (i = 0; i < 256; ++i) cnt[i] = 0;
    for (k = 0; k < (1U<<(ht)->bits); ++k)
        if (ht->values[k] != 0)
            ++cnt[ht->values[k] < 256 ? ht->values[k] : 255];
    for (i = 1; i < 256; ++i)
        printf("%d\t%ld\n", i, (long)cnt[i]);	
}

static int count_file(const char *fn, int k, unsigned int p)
{
	//gzFile fp;
	//kseq_t *ks;
	HashTable *ht;
	unsigned int capacity = 1U << p;
    unsigned int read_count = 0;



    // variables para cuda
	//HashTable *ht_d;
	//char **reads_d;
	//unsigned int *cnt_d;

    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(fn, "r");
    if (fp == NULL) exit(EXIT_FAILURE);

   
	ReadSeqList *current, *head;
	head = current = NULL;

    while ((read = getline(&line, &len, fp)) != -1) {

        read_count++;

		ReadSeqList *node = (ReadSeqList*)malloc(sizeof(ReadSeqList));
        node->sequence = (char*)malloc(strlen(line));
        strcpy(node->sequence, line);
        node->length = read;
        node->next =NULL;

        if(head == NULL){
            current = head = node;
        } else {
            current = current->next = node;
        }


    }

    fclose(fp);
    if (line) free(line);

    // leer los reads
	//if ((fp = gzopen(fn, "r")) == 0) return 0;
	//ks = kseq_init(fp); // descriptor fichero fastaq

	// leer los datos del fichero de entrada y guardarlos en memoria
/*
	while (kseq_read(ks) >= 0) {
        read_count++;

		ReadSeqList *node = (ReadSeqList*)malloc(sizeof(ReadSeqList));
        node->sequence = (char*)malloc(strlen(ks->seq.s) + 1);
        strcpy(node->sequence, ks->seq.s);
        node->length = ks->seq.l;
        node->next =NULL;

        if(head == NULL){
            current = head = node;
        } else {
            current = current->next = node;
        }

	}
*/
	//kseq_destroy(ks);
	//gzclose(fp);

    // crear un array de tamaño fijo para almacenar las lecturas
    // char **reads = malloc(read_count * sizeof(char*));



    // inicializar hashtable
	ht = HashTable_init(p, read_count);





	// allocate memory in device
	//cudaMalloc((void **)&reads_d, read_count * sizeof(char *));
	//cudaMalloc((void **)&ht_d, sizeof(HashTable));
	//cudaMalloc((void **)ht_d->keys, capacity * sizeof(unsigned long long int));
	//cudaMalloc((void **)ht_d->values, capacity * sizeof(unsigned int));
	//cudaMalloc((void **)ht_d->collition, capacity * sizeof(unsigned int));
	//cudaMalloc((void **)&cnt_d, 256 * sizeof(unsigned int));

	//cudaMemset(ht_d->keys, 0, capacity * sizeof(unsigned long long int));
	//cudaMemset(ht_d->values, 0, capacity * sizeof(unsigned int));
	//cudaMemset(cnt_d, 0, 256 * sizeof(unsigned int));

	// copy data to device
	//cudaMemcpy(ht_d, ht, sizeof(HashTable), cudaMemcpyHostToDevice);





	// copiar los read a la memoria de la GPU
	//char **temp_reads_d = (char **)malloc(read_count * sizeof(char *));
    //unsigned int i;
    char** reads = malloc(read_count * sizeof(char*));
    unsigned int i;
	for(i=0, current = head; current; current=current->next){
        reads[i] = current->sequence;
		i++;
    }



    //printf("%s\n", reads[0]);

    printf("total reads: %d\n", read_count);


    // invocar kernels
    //int thr = 1024;

    current = head;
    for(i = 0; i<read_count ; i++){
        kernel_count_seq_kmers(ht, k, reads[i]);
    }

	kernel_print_hist(ht);



	//cudaMemcpy(ht, ht_d, sizeof(HashTable), cudaMemcpyDeviceToHost);
	//cudaMemcpy(ht->keys, ht_d->keys, capacity * sizeof(unsigned long long int), cudaMemcpyDeviceToHost);
	//cudaMemcpy(ht->values, ht_d->values, capacity * sizeof(unsigned int), cudaMemcpyDeviceToHost);
	//cudaMemcpy(cnt, cnt_d, capacity * sizeof(unsigned int), cudaMemcpyDeviceToHost);


    printf("COUNT: %d\n\n", ht->count);

    // for(i = 0; i< 5000 ; i++){
    //     printf("key: ");
    //     printf("%"PRIu64"\t", ht->keys[i]);
    //     printf("value: %d\t", ht->values[i]);
    //     printf("collitions: %d\n", ht->collition[i]);
    // }


	// limpieza
	for(current = head; current; current=current->next){
        free(current->sequence);
        free(current);
    }

	//cudaFree(reads_d);
	//cudaFree(ht_d);
	//cudaFree(cnt_d);

	HashTable_destory(ht);

    return 0;
}


int main(int argc, char *argv[])
{
	int k = 31;
    unsigned int p = 27;

/*
	ketopt_t o = KETOPT_INIT;
	while ((c = ketopt(&o, argc, argv, 1, "k:", 0)) >= 0)
		if (c == 'k') k = atoi(o.arg);
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: kc-c1 [-k %d] <in.fa>\n", k);
		return 1;
	}
*/
    k = (int)strtol(argv[1], NULL, 10);
    p = (unsigned int)strtol(argv[2], NULL, 10);
	count_file(argv[3], k, p);

	return 0;
}
