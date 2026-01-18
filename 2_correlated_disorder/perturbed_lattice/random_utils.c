#include "random_utils.h"

// Xoshiro256+ implementation for better RNG
void init_rng(RNG* rng, unsigned long seed) {
    // Initialize all state variables from seed for reproducibility
    rng->state[0] = seed ^ 0x123456789ABCDEF0UL;
    rng->state[1] = (seed << 1) ^ 0xFEDCBA9876543210UL;
    rng->state[2] = (seed << 2) ^ 0x0123456789ABCDEFUL;
    rng->state[3] = (seed << 3) ^ 0xDEADBEEFC0FFEE00UL;
    
    // Mix the state to avoid poor initial values
    for (int i = 0; i < 20; i++) {
        uniform_random(rng);
    }
}

static inline unsigned long rotl(const unsigned long x, int k) {
    return (x << k) | (x >> (64 - k));
}

double uniform_random(RNG* rng) {
    const unsigned long result = rng->state[0] + rng->state[3];
    const unsigned long t = rng->state[1] << 17;

    rng->state[2] ^= rng->state[0];
    rng->state[3] ^= rng->state[1];
    rng->state[1] ^= rng->state[2];
    rng->state[0] ^= rng->state[3];

    rng->state[2] ^= t;
    rng->state[3] = rotl(rng->state[3], 45);

    return (result >> 11) * 0x1.0p-53;
}

double gaussian_random(RNG* rng, double mu, double sigma) {
    static int has_spare = 0;
    static double spare;
    
    if (has_spare) {
        has_spare = 0;
        return spare * sigma + mu;
    }
    
    has_spare = 1;
    double u, v, s;
    do {
        u = uniform_random(rng) * 2.0 - 1.0;
        v = uniform_random(rng) * 2.0 - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    
    double mult = sqrt(-2.0 * log(s) / s);
    spare = v * mult;
    return mu + sigma * u * mult;
}
