#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fq.h>
#include <flint/fmpz_poly.h>
#include <flint/ulong_extras.h>
#include <flint/fq_poly.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>
#include <sys/time.h>

const bool DEBUG = false;
const bool TIMER = true;
const bool PRINT_CFF = false;
const bool VERIFY_CFF = false;

double get_elapsed_ms(struct timeval start_time, struct timeval end_time) {
    return (end_time.tv_sec - start_time.tv_sec) * 1000.0 + (end_time.tv_usec - start_time.tv_usec) / 1000.0;
}

static inline ulong get_matrix_pos(ulong i, ulong j, ulong num_rows) {
    return j * num_rows + i;
}

void initialize_finite_field(ulong p, slong n, fq_ctx_t ctx, fmpz_t p_fmpz) {
    fmpz_init_set_ui(p_fmpz, p);
    fq_ctx_init(ctx, p_fmpz, n, "x");
}

void create_finite_field_elements(ulong p, ulong n, ulong q, fq_t field_elements[], fq_ctx_t ctx) {
    struct timeval start_time, end_time;
    if (TIMER) gettimeofday(&start_time, NULL);
    if (DEBUG) flint_printf("size of finite field: %u\n", q);
    fmpz_poly_t poly;
    for (ulong i = 0; i < q; i++) {
        // Creates the polynomial representation of the element
        fmpz_poly_init(poly);
        ulong temp = i;
        for (ulong j = 0; j < n; j++) {
            fmpz_poly_set_coeff_ui(poly, j, temp % p);
            temp /= p;
        }

        fq_init(field_elements[i], ctx);
        fq_set_fmpz_poly(field_elements[i], poly, ctx);

        fmpz_poly_clear(poly);

        if (DEBUG) {
            flint_printf("\nelement index: %u\n", i);
            printf("element pretty print: ");
            fq_print_pretty(field_elements[i], ctx);
            printf("\n");
            printf("element raw print: ");
            fq_print(field_elements[i], ctx);
            printf("\n");
        }
    }
    if (DEBUG) printf("\n\n");
    if (TIMER) {
        gettimeofday(&end_time, NULL);
        printf("Elapsed time to create finite field elements: %f ms\n", get_elapsed_ms(start_time, end_time));
    }
}

bool* create_cff_matrix(ulong q, slong k, fq_t field_elements[], fq_ctx_t ctx) {
    struct timeval start_time, end_time;
    if (TIMER) gettimeofday(&start_time, NULL);

    ulong num_rows = q * q;
    ulong num_polys = n_pow(q,k+1);

    bool* cff = (bool*)malloc(num_rows * num_polys * sizeof(bool));
    if (cff == NULL) {
        flint_printf("Failed to allocate memory for CFF matrix with dimensions %lu x %lu\n", q * q, num_polys);
        exit(EXIT_FAILURE);
    }

    #pragma omp parallel for shared(cff, field_elements, q, num_polys, ctx)
    for (ulong j = 0; j < num_polys; j++) {
        // Creates the polynomial 
        fq_poly_t polynomial;
        fq_poly_init(polynomial, ctx);
        ulong temp = j;
        for (ulong m = 0; m <= k; m++) {
            fq_poly_set_coeff(polynomial, m, field_elements[temp % q], ctx);
            temp /= q;
        }

        // Evaluates the polynomial in all elements of the finite field
        fq_t fx;
        ulong x_i = 0;
        ulong fx_i = 0;
        for (ulong i = 0; i < q * q; i += q) {
            fq_init(fx, ctx);
            fq_poly_evaluate_fq(fx, polynomial, field_elements[x_i], ctx);
            for (ulong k = 0; k < q; k++) {
                cff[get_matrix_pos(i + k, j, num_rows)] = fq_equal(fx, field_elements[fx_i], ctx);
                fx_i++;
            }
            x_i++;
            fx_i = 0;
            fq_clear(fx, ctx);
        }

        fq_poly_clear(polynomial, ctx);
    }

    if (TIMER) {
        gettimeofday(&end_time, NULL);
        printf("Elapsed time to create CFF matrix: %f ms\n", get_elapsed_ms(start_time, end_time));
    }
    return cff;
}


void clear_memory(ulong q, fq_t field_elements[], bool* cff, fq_ctx_t ctx, fmpz_t p_fmpz) {
    for (ulong i = 0; i < q; i++) {
        fq_clear(field_elements[i], ctx);
    }
    fq_ctx_clear(ctx);
    fmpz_clear(p_fmpz);
}

bool* get_cff(ulong p, slong n, slong k) {
    struct timeval start_time, end_time;
    if (TIMER) gettimeofday(&start_time, NULL);
    fq_ctx_t ctx;
    fmpz_t p_fmpz;
    initialize_finite_field(p, n, ctx, p_fmpz);

    ulong q = n_pow(p, n);
    fq_t field_elements[q];
    create_finite_field_elements(p, n, q, field_elements, ctx);
    
    bool* cff = create_cff_matrix(q, k, field_elements, ctx);
    
    clear_memory(q, field_elements, cff, ctx, p_fmpz);

    if (TIMER) {
        gettimeofday(&end_time, NULL);
        printf("Total elapsed time: %f ms\n", get_elapsed_ms(start_time, end_time));
    }
    return cff;
}

void print_cff(bool* cff, ulong q, ulong k) {
    ulong num_polys = n_pow(q,k+1);
    ulong num_rows = q * q;
    printf("CFF:\n");
    for (ulong i = 0; i < q * q; i++) {
        for (ulong j = 0; j < num_polys; j++) {
            printf("%d ", cff[get_matrix_pos(i, j, num_rows)]);
        }
        printf("\n");
    }
}

bool have_unique_element(bool* cff, ulong num_lines, ulong num_collumns, ulong current_columns[], ulong comb_size) {
    bool found[comb_size];
    for (ulong k = 0; k < comb_size; k++) {
        found[k] = false;
    }
    for (ulong i = 0; i < num_lines; i++) {
        ulong found_k = -1;
        for (ulong k = 0; k < comb_size; k++) {
            ulong j = current_columns[k];
            if (cff[get_matrix_pos(i, j, num_lines)] == 1) {
                if (found_k == -1) {
                    found_k = k;
                } else {
                    found_k = -1;
                    break;
                }
            } 
        }
        if (found_k != -1) {
            found[found_k] = true;
        }
    }
    for (ulong k = 0; k < comb_size; k++) {
        if (!found[k]) return false;
    }
    return true;
}

void print_combination(ulong current_columns[], ulong comb_size) {
    for (ulong i = 0; i < comb_size; i++) {
        printf("%lu ", current_columns[i]);
    }
    printf("\n");
}

bool is_cff(ulong q, ulong k, ulong d, bool* cff) {
    ulong num_collumns = n_pow(q,k+1);
    ulong num_lines = q * q;
    ulong comb_size = d + 1;
    if (comb_size > num_collumns) return false;

    ulong current_columns[comb_size];
    for (ulong i = 0; i < comb_size; i++) {
        current_columns[i] = i;
    }

    while (true) {
        if (!have_unique_element(cff, num_lines, num_collumns, current_columns, comb_size)) {
            printf("Failed columns combination: ");
            print_combination(current_columns, comb_size);
            return false;
        }

        ulong next_increment_index = comb_size - 1;
        ulong max_value = num_collumns - (comb_size - next_increment_index);
        while (next_increment_index >= 0 && current_columns[next_increment_index] == max_value) {
            next_increment_index--;
            max_value = num_collumns - (comb_size - next_increment_index);
        }

        if (next_increment_index < 0) {
            break;
        }

        current_columns[next_increment_index]++;

        // Reset all subsequent indices
        for (ulong i = next_increment_index + 1; i < comb_size; i++) {
            current_columns[i] = current_columns[i - 1] + 1;
        }

    }
    return true;
}

int main() {
    const ulong p = 2;
    const slong n = 3;
    const slong k = 4;

    ulong q = n_pow(p, n);

    bool* cff = get_cff(p, n, k);
    printf("CFF generated\n");

    if (PRINT_CFF) print_cff(cff, q, k);

    if (VERIFY_CFF) {
        const ulong d = 1;
        if (is_cff(q, k, d, cff)) {
            flint_printf("It is a %u-CFF\n", d);
        } else {
            flint_printf("It is not a %u-CFF\n", d);
        }
    }

    free(cff);

    return 0;
}