/*
Benchmark code as suggested in "How to Benchmark Code Execution Times on 
Intel® IA - 32 and IA - 64 Instruction Set Architectures"

Handles serialization and out-of-order execution

arguments are 32-bit unsigned ints (uint32)

For cleaner usage, the START and STOP macros assume that following are declared:
unsigned rdtsc_start_high, rdtsc_start_low, rdtsc_stop_high,rdtsc_top_low
uint64_t rdtsc_count[N] , with with N chosen to provide sufficient storage
unsigned rdtsc_idx

Then, usage becomes:
RDTSC_START
..do stuff..
RDTSC_STOP

print rdtsc_count[i] to see TSC counts elapsed

RDTSC_RAW only does a RDTSC, without queue flushes.  Can be used
for checkpointing with minimal overhead
*/

#define RDTSC_START \
	__asm volatile ( \
		"CPUID\n\t" \
		"RDTSC\n\t" \
		"mov %%edx, %0\n\t" \
		"mov %%eax, %1\n\t": "=r" (rdtsc_start_high), "=r" (rdtsc_start_low):: \
			"%rax", "%rbx", "%rcx", "%rdx");

#define RDTSC_STOP \
	__asm volatile ( \
		"RDTSCP\n\t" \
		"mov %%edx, %0\n\t" \
		"mov %%eax, %1\n\t" \
		"CPUID\n\t" : "=r" (rdtsc_stop_high), "=r" (rdtsc_stop_low):: \
			"%rax", "%rbx", "%rcx", "%rdx"); \
	rdtsc_counts[rdtsc_idx] = (((uint64_t)rdtsc_stop_high  << 32) | rdtsc_stop_low  ) - \
	                    (((uint64_t)rdtsc_start_high << 32) | rdtsc_start_low ); \
	++rdtsc_idx;

#define RDTSC_RAW(low,high) \
	__asm volatile ( \
		"RDTSC\n\t" \
		"mov %%edx, %0\n\t" \
		"mov %%eax, %1\n\t": "=r" (high), "=r" (low):: \
			"%eax", "%edx");
