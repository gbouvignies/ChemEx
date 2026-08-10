#include <immintrin.h>

unsigned int chemex_get_mxcsr(void) {
    return _mm_getcsr();
}

unsigned short chemex_get_x87_control_word(void) {
    unsigned short control_word;
    __asm__ volatile("fnstcw %0" : "=m"(control_word));
    return control_word;
}
