include ../../config_defs

all: pmsort_test pssort_test pssort_bench pssort_bench_record

CFLAGS = -O3 -g
# -qarch=450d -qtune=450 

# position-independent code
ifeq ($(FPIC), YES)
CFLAGS += -fPIC
endif

pssort.o: pssort.c
	$(CC) $(CFLAGS) -c pssort.c

pmsort.o: pmsort.c
	$(CC) $(CFLAGS) -c pmsort.c

util.o: util.c
	$(CC) $(CFLAGS) -c util.c

pssort_test.o: pssort_test.c
	$(CC) $(CFLAGS) -c pssort_test.c

pssort_test: util.o pssort.o pmsort.o pssort_test.o
	$(CC) $(CFLAGS) -o pssort_test util.o pssort.o pmsort.o pssort_test.o

pmsort_test.o: pmsort_test.c
	$(CC) $(CFLAGS) -c pmsort_test.c

pmsort_test: util.o pmsort.o pmsort_test.o
	$(CC) $(CFLAGS) -o pmsort_test util.o pmsort.o pmsort_test.o

pssort_bench.o: pssort_bench.c
	$(CC) $(CFLAGS) -c pssort_bench.c

pssort_bench: util.o pssort.o pmsort.o pssort_bench.o
	$(CC) $(CFLAGS) -o pssort_bench util.o pssort.o pmsort.o pssort_bench.o

pssort_bench_record.o: pssort_bench_record.c
	$(CC) $(CFLAGS) -c pssort_bench_record.c

pssort_bench_record: util.o pssort.o pmsort.o pssort_bench_record.o
	$(CC) $(CFLAGS) -o pssort_bench_record util.o pssort.o pmsort.o pssort_bench_record.o

clean:
	rm -f *.o pssort_test pmsort_test pssort_bench pssort_bench_record


