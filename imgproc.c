#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "proto.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <dirent.h>    
#include <time.h>  
#include <arm_neon.h>

double total_time_program = 0.0;
double total_time_laplacian = 0.0;
Image *Laplacian(Image *);
Image *Sobel(Image *);
Image *ReadPNMImage(char *);
Image *CreateNewImage(Image *, char *comment);
int boundaryCheck(int index_x, int index_y, int width, int height);
int TestReadImage(char *, char *);
void SavePNMImage(Image *, char *, char*, char*);
char *Extract_Filename_Stem(char*);
void process_folder_to_pgm_then_run(char *input_dir, char *pgm_output_dir);
double now_seconds(void);

int main(int argc, char **argv)
{
    char *input_dir = (argc > 1) ? argv[1] : "images";
    char *pgm_output_dir = (argc > 2) ? argv[2] : "grayscale_inputs_pgm";
    double start = now_seconds();
    process_folder_to_pgm_then_run(input_dir, pgm_output_dir);
    double end = now_seconds();
    total_time_program += end - start;
    
    printf("\n###############################\n");
    printf("#                             #\n");
    printf("#                             #\n");
    printf("### Program time: %.3f s ####\n", total_time_program);
    printf("# Laplacian operator time: %.3f s #\n", total_time_laplacian);
    printf("#                             #\n");
    printf("#                             #\n");
    printf("###############################\n\n");

    return(0);
}
double now_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);                 // C11 standard
    return ts.tv_sec + ts.tv_nsec / 1e9;
}

void ensure_output_dir(const char *dir) {
    struct stat st;

    if (stat(dir, &st) == -1) {
        if (mkdir(dir, 0755) == -1) {
            perror("mkdir failed");
            exit(1);
        }
    }
}

static int convert_to_pgm(char *input_path, char *output_dir, char *out_path, size_t out_path_sz) {
    ensure_output_dir(output_dir);

    char *stem = Extract_Filename_Stem(input_path);
    if (!stem) {
        fprintf(stderr, "OOM extracting stem for %s\n", input_path);
        return 1;
    }

    snprintf(out_path, out_path_sz, "%s/%s.pgm", output_dir, stem);

    char cmd[2048];
    snprintf(cmd, sizeof(cmd), "magick \"%s\" -colorspace Gray -depth 8 \"%s\"", input_path, out_path);

    int rc = system(cmd);
    if (rc != 0) {
        fprintf(stderr, "ImageMagick conversion failed (%d): %s\n", rc, input_path);
        free(stem);
        return rc;
    }

    printf("Converted -> %s\n", out_path);
    free(stem);
    return 0;
}

// Traverse input_dir, convert each supported image to PGM
void process_folder_to_pgm_then_run(char *input_dir, char *pgm_output_dir) {
    DIR *dir = opendir(input_dir);
    if (!dir) {
        perror("opendir failed");
        exit(1);
    }

    struct dirent *ent;
    while ((ent = readdir(dir)) != NULL) {
        if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0) continue;
        char in_path[1024];
        snprintf(in_path, sizeof(in_path), "%s/%s", input_dir, ent->d_name);
        char pgm_path[1024];
        if (convert_to_pgm(in_path, pgm_output_dir, pgm_path, sizeof(pgm_path)) == 0) {
            // Now run your pipeline on the generated PGM
            TestReadImage(pgm_path, "");
        }
    }

    closedir(dir);
}

int TestReadImage(char *file_in, char *file_out)
{
    Image *image;
    Image *laplacian, *sobel;
    char *stem = Extract_Filename_Stem(file_in);
    image = ReadPNMImage(file_in);
    double start = now_seconds();
    laplacian = Laplacian(image);
    // sobel = Sobel(image);
    double end = now_seconds();
    total_time_laplacian += (end - start);
    char filename[300];
    snprintf(filename, sizeof(filename), "%s.pgm", stem);

    SavePNMImage(laplacian, filename, "laplacian_pgm", "laplacian_png");
    // SavePNMImage(sobel, filename, "sobel_pgm", "sobel_png");
    free(stem);
    return(0);
}

static inline unsigned char clamp_u8_int(int v) {
    if (v < 0)   return 0;
    if (v > 255) return 255;
    return (unsigned char)v;
}

static inline uint8x16_t load_u8_zeropad_right(const uint8_t *p, int n_rem) {
    uint8_t tmp[16] = {0};
    if (n_rem > 0) memcpy(tmp, p, (size_t)n_rem);
    return vld1q_u8(tmp);
}

Image *Laplacian(Image *image) {
    const int W = image->Width;
    const int H = image->Height;

    unsigned char *src = image->data;
    Image *out = CreateNewImage(image, (char*)"#testing function");
    unsigned char *dst = out->data;

// --- TOP ROW (y = 0): sum = bottom_left + bottom_right ---
if (H > 0) {
    const int y = 0;
    const uint8_t *rowC = src + y*W;
    const uint8_t *rowB = (H > 1) ? (src + (y+1)*W) : NULL; // may be NULL if H==1
    uint8_t *rowO = dst + y*W;

    if (W >= 16 && rowB) {
        int x = 0;

        // prevB starts as zeros to handle left edge for left-shift-by-1
        uint8x16_t prevB = vdupq_n_u8(0);

        // process full 16-byte blocks; we also need the next block to build right-shift-by-1
        for (; x + 16 <= W; x += 16) {
            // current and next (for right shift); next is zero-padded at the very end
            uint8x16_t curC = vld1q_u8(rowC + x);
            uint8x16_t curB = vld1q_u8(rowB + x);
            uint8x16_t nextB;
            if (x + 32 <= W) {
                nextB = vld1q_u8(rowB + x + 16);
            } else {
                int rem = W - (x + 16);
                nextB = load_u8_zeropad_right(rowB + x + 16, rem > 0 ? rem : 0);
            }

            // left_b = rowB shifted right by 1: [x-1 .. x+14], pad 0 for x=0 (handled by prevB)
            uint8x16_t left_b  = vextq_u8(prevB, curB, 15);

            // right_b = rowB shifted left by 1:  [x+1 .. x+16], pad 0 after last pixel (via nextB)
            uint8x16_t right_b = vextq_u8(curB, nextB, 1);

            // widen to u16 and sum
            uint8x8_t  lb_lo = vget_low_u8(left_b),  lb_hi = vget_high_u8(left_b);
            uint8x8_t  rb_lo = vget_low_u8(right_b), rb_hi = vget_high_u8(right_b);
            uint16x8_t sum_lo = vaddl_u8(lb_lo, rb_lo);
            uint16x8_t sum_hi = vaddl_u8(lb_hi, rb_hi);

            // 4 * center
            uint8x8_t  c_lo = vget_low_u8(curC), c_hi = vget_high_u8(curC);
            uint16x8_t four_c_lo = vshll_n_u8(c_lo, 2);
            uint16x8_t four_c_hi = vshll_n_u8(c_hi, 2);

            // res = 4c - sum (signed 16) -> clamp -> u8
            int16x8_t res_lo = vsubq_s16(vreinterpretq_s16_u16(four_c_lo),
                                         vreinterpretq_s16_u16(sum_lo));
            int16x8_t res_hi = vsubq_s16(vreinterpretq_s16_u16(four_c_hi),
                                         vreinterpretq_s16_u16(sum_hi));
            uint8x8_t out_lo = vqmovun_s16(res_lo);
            uint8x8_t out_hi = vqmovun_s16(res_hi);
            uint8x16_t out_px = vcombine_u8(out_lo, out_hi);

            vst1q_u8(rowO + x, out_px);

            // advance prevB for the next block
            prevB = curB;
        }

        // tail (scalar) for any leftover pixels
        for (; x < W; ++x) {
            int sum = 0;
            if (x - 1 >= 0) sum += rowB[x - 1];
            if (x + 1 <  W) sum += rowB[x + 1];
            int v = 4*rowC[x] - sum;
            rowO[x] = clamp_u8_int(v);
        }
    } else {
        // scalar fallback (W < 16 or no rowB)
        for (int x = 0; x < W; ++x) {
            int sum = 0;
            if (rowB) {
                if (x - 1 >= 0) sum += rowB[x - 1];
                if (x + 1 <  W) sum += rowB[x + 1];
            }
            int v = 4*rowC[x] - sum;
            rowO[x] = clamp_u8_int(v);
        }
    }
}

if (H > 1) {
    const int y = H - 1;
    const uint8_t *rowC = src + y*W;
    const uint8_t *rowT = src + (y-1)*W;
    uint8_t *rowO = dst + y*W;

    if (W >= 16) {
        int x = 0;
        uint8x16_t prevT = vdupq_n_u8(0);

        for (; x + 16 <= W; x += 16) {
            uint8x16_t curC = vld1q_u8(rowC + x);
            uint8x16_t curT = vld1q_u8(rowT + x);
            uint8x16_t nextT;
            if (x + 32 <= W) {
                nextT = vld1q_u8(rowT + x + 16);
            } else {
                int rem = W - (x + 16);
                nextT = load_u8_zeropad_right(rowT + x + 16, rem > 0 ? rem : 0);
            }

            // top-left  = rowT shifted right by 1
            // top-right = rowT shifted left  by 1
            uint8x16_t tl = vextq_u8(prevT, curT, 15);
            uint8x16_t tr = vextq_u8(curT,  nextT, 1);

            uint8x8_t  tl_lo = vget_low_u8(tl), tl_hi = vget_high_u8(tl);
            uint8x8_t  tr_lo = vget_low_u8(tr), tr_hi = vget_high_u8(tr);
            uint16x8_t sum_lo = vaddl_u8(tl_lo, tr_lo);
            uint16x8_t sum_hi = vaddl_u8(tl_hi, tr_hi);

            uint8x8_t  c_lo = vget_low_u8(curC), c_hi = vget_high_u8(curC);
            uint16x8_t four_c_lo = vshll_n_u8(c_lo, 2);
            uint16x8_t four_c_hi = vshll_n_u8(c_hi, 2);

            int16x8_t res_lo = vsubq_s16(vreinterpretq_s16_u16(four_c_lo),
                                         vreinterpretq_s16_u16(sum_lo));
            int16x8_t res_hi = vsubq_s16(vreinterpretq_s16_u16(four_c_hi),
                                         vreinterpretq_s16_u16(sum_hi));
            uint8x8_t out_lo = vqmovun_s16(res_lo);
            uint8x8_t out_hi = vqmovun_s16(res_hi);
            uint8x16_t out_px = vcombine_u8(out_lo, out_hi);

            vst1q_u8(rowO + x, out_px);

            prevT = curT;
        }

        for (; x < W; ++x) {
            int sum = 0;
            if (x - 1 >= 0) sum += rowT[x - 1];
            if (x + 1 <  W) sum += rowT[x + 1];
            int v = 4*rowC[x] - sum;
            rowO[x] = clamp_u8_int(v);
        }
    } else {
        // scalar fallback for narrow rows
        for (int x = 0; x < W; ++x) {
            int sum = 0;
            if (x - 1 >= 0) sum += rowT[x - 1];
            if (x + 1 <  W) sum += rowT[x + 1];
            int v = 4*rowC[x] - sum;
            rowO[x] = clamp_u8_int(v);
        }
    }
}

    // --- SIMD interior rows: y in [1, H-2] ---
    for (int y = 1; y <= H - 2; ++y) {
        const unsigned char *rowT = src + (y-1)*W;
        const unsigned char *rowC = src +  y   *W;
        const unsigned char *rowB = src + (y+1)*W;
        unsigned char       *rowO = dst +  y   *W;

        // Left border (x=0) scalar
        {
            int x = 0;
            int sum = 0;
            if (x-1 >= 0) { sum += rowT[x-1]; sum += rowB[x-1]; }
            if (x+1 <  W) { sum += rowT[x+1]; sum += rowB[x+1]; }
            int v = 4*rowC[x] - sum;
            rowO[x] = clamp_u8_int(v);
        }

        // SIMD for safe interior: ensure we can read (x-1) .. (x+16) and (x+1) .. (x+17)
        // That means x runs up to W-18 inclusive (so x+16 <= W-2 and x+1+16 <= W-1)
        int x = 1;
        int simd_limit = (W >= 18) ? (W - 18) : 1;
        for (; x <= simd_limit; x += 16) {
            // Load centers
            uint8x16_t c  = vld1q_u8(rowC + x);

            // Load diagonals (unaligned ok on ARM64):
            // top-left  at x-1 .. x+14
            // top-right at x+1 .. x+16
            // bot-left  at x-1 .. x+14
            // bot-right at x+1 .. x+16
            uint8x16_t tl = vld1q_u8(rowT + (x - 1));
            uint8x16_t tr = vld1q_u8(rowT + (x + 1));
            uint8x16_t bl = vld1q_u8(rowB + (x - 1));
            uint8x16_t br = vld1q_u8(rowB + (x + 1));

            // Widen to 16-bit
            uint8x8_t  c_lo  = vget_low_u8(c),   c_hi  = vget_high_u8(c);
            uint8x8_t  tl_lo = vget_low_u8(tl),  tl_hi = vget_high_u8(tl);
            uint8x8_t  tr_lo = vget_low_u8(tr),  tr_hi = vget_high_u8(tr);
            uint8x8_t  bl_lo = vget_low_u8(bl),  bl_hi = vget_high_u8(bl);
            uint8x8_t  br_lo = vget_low_u8(br),  br_hi = vget_high_u8(br);

            uint16x8_t c16_lo  = vmovl_u8(c_lo);
            uint16x8_t c16_hi  = vmovl_u8(c_hi);
            uint16x8_t tl16_lo = vmovl_u8(tl_lo), tl16_hi = vmovl_u8(tl_hi);
            uint16x8_t tr16_lo = vmovl_u8(tr_lo), tr16_hi = vmovl_u8(tr_hi);
            uint16x8_t bl16_lo = vmovl_u8(bl_lo), bl16_hi = vmovl_u8(bl_hi);
            uint16x8_t br16_lo = vmovl_u8(br_lo), br16_hi = vmovl_u8(br_hi);

            // sum of diagonals
            uint16x8_t sum_lo = vaddq_u16(vaddq_u16(tl16_lo, tr16_lo), vaddq_u16(bl16_lo, br16_lo));
            uint16x8_t sum_hi = vaddq_u16(vaddq_u16(tl16_hi, tr16_hi), vaddq_u16(bl16_hi, br16_hi));

            // 4 * center
            uint16x8_t four_c_lo = vshlq_n_u16(c16_lo, 2);
            uint16x8_t four_c_hi = vshlq_n_u16(c16_hi, 2);

            // res = 4c - sum (in signed 16)
            int16x8_t res_lo = vsubq_s16(vreinterpretq_s16_u16(four_c_lo),
                                         vreinterpretq_s16_u16(sum_lo));
            int16x8_t res_hi = vsubq_s16(vreinterpretq_s16_u16(four_c_hi),
                                         vreinterpretq_s16_u16(sum_hi));

            // Saturating narrow to u8 (clamp to [0,255])
            uint8x8_t out_lo = vqmovun_s16(res_lo);
            uint8x8_t out_hi = vqmovun_s16(res_hi);
            uint8x16_t out_px = vcombine_u8(out_lo, out_hi);

            vst1q_u8(rowO + x, out_px);
        }

        // Remaining interior pixels up to W-2 (scalar)
        for (; x <= W - 2; ++x) {
            int sum = rowT[x-1] + rowT[x+1] + rowB[x-1] + rowB[x+1];
            int v = 4*rowC[x] - sum;
            rowO[x] = clamp_u8_int(v);
        }

        // Right border (x=W-1) scalar
        if (W >= 2) {
            int xb = W - 1;
            int sum = 0;
            sum += rowT[xb-1];            // (y-1, x-1)
            // (y-1, x+1) out of bounds -> +0
            sum += rowB[xb-1];            // (y+1, x-1)
            // (y+1, x+1) out of bounds -> +0
            int v = 4*rowC[xb] - sum;
            rowO[xb] = clamp_u8_int(v);
        }
    }

    return out;
}

Image *Sobel(Image *image) {
    const int W = image->Width;
    const int H = image->Height;

    const uint8_t *src = image->data;
    Image *out = CreateNewImage(image, (char*)"#Sobel NEON");
    uint8_t *dst = out->data;
    
    for (int y = 0; y < H; ++y) {
        const uint8_t *rowT = (y > 0)     ? (src + (y-1)*W) : NULL;
        const uint8_t *rowC =               src +  y   *W;
        const uint8_t *rowB = (y+1 < H)   ? (src + (y+1)*W) : NULL;
        uint8_t       *rowO =               dst +  y   *W;

        if (W >= 16) {
            int x = 0;

            // For vext shifts we need previous & next 16B blocks per row.
            uint8x16_t prevT = vdupq_n_u8(0), prevC = vdupq_n_u8(0), prevB = vdupq_n_u8(0);

            for (; x + 16 <= W; x += 16) {
                // current blocks
                uint8x16_t curT = rowT ? vld1q_u8(rowT + x) : vdupq_n_u8(0);
                uint8x16_t curC =        vld1q_u8(rowC + x);
                uint8x16_t curB = rowB ? vld1q_u8(rowB + x) : vdupq_n_u8(0);

                // next blocks (for right shift), zero-padded at image tail
                uint8x16_t nextT, nextC, nextB;
                if (x + 32 <= W) {
                    nextT = rowT ? vld1q_u8(rowT + x + 16) : vdupq_n_u8(0);
                    nextC =        vld1q_u8(rowC + x + 16);
                    nextB = rowB ? vld1q_u8(rowB + x + 16) : vdupq_n_u8(0);
                } else {
                    int rem = W - (x + 16);
                    nextT = rowT ? load_u8_zeropad_right(rowT + x + 16, rem > 0 ? rem : 0) : vdupq_n_u8(0);
                    nextC =        load_u8_zeropad_right(rowC + x + 16, rem > 0 ? rem : 0);
                    nextB = rowB ? load_u8_zeropad_right(rowB + x + 16, rem > 0 ? rem : 0) : vdupq_n_u8(0);
                }

                // Left/Right neighbors via shifts
                // left = shift-right by 1  (uses prev as carry-in)
                // right= shift-left  by 1  (uses next as carry-in)
                uint8x16_t lT = vextq_u8(prevT, curT, 15);
                uint8x16_t rT = vextq_u8(curT,  nextT, 1);
                uint8x16_t lC = vextq_u8(prevC, curC, 15);
                uint8x16_t rC = vextq_u8(curC,  nextC, 1);
                uint8x16_t lB = vextq_u8(prevB, curB, 15);
                uint8x16_t rB = vextq_u8(curB,  nextB, 1);

                // Widen to s16 (signed) for arithmetic
                int16x8_t lT_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(lT)));
                int16x8_t lT_hi = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(lT)));
                int16x8_t rT_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(rT)));
                int16x8_t rT_hi = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(rT)));

                int16x8_t lC_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(lC)));
                int16x8_t lC_hi = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(lC)));
                int16x8_t rC_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(rC)));
                int16x8_t rC_hi = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(rC)));

                int16x8_t lB_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(lB)));
                int16x8_t lB_hi = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(lB)));
                int16x8_t rB_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(rB)));
                int16x8_t rB_hi = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(rB)));

                int16x8_t cT_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(curT)));
                int16x8_t cT_hi = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(curT)));
                int16x8_t cC_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(curC)));
                int16x8_t cC_hi = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(curC)));
                int16x8_t cB_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(curB)));
                int16x8_t cB_hi = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(curB)));

                // Gx = (rT + 2*rC + rB) - (lT + 2*lC + lB)
                int16x8_t rx_lo = vaddq_s16(vaddq_s16(rT_lo, vshlq_n_s16(rC_lo, 1)), rB_lo);
                int16x8_t rx_hi = vaddq_s16(vaddq_s16(rT_hi, vshlq_n_s16(rC_hi, 1)), rB_hi);
                int16x8_t lx_lo = vaddq_s16(vaddq_s16(lT_lo, vshlq_n_s16(lC_lo, 1)), lB_lo);
                int16x8_t lx_hi = vaddq_s16(vaddq_s16(lT_hi, vshlq_n_s16(lC_hi, 1)), lB_hi);
                int16x8_t gx_lo = vsubq_s16(rx_lo, lx_lo);
                int16x8_t gx_hi = vsubq_s16(rx_hi, lx_hi);

                // Gy = (lB + 2*cB + rB) - (lT + 2*cT + rT)
                int16x8_t by_lo = vaddq_s16(vaddq_s16(lB_lo, vshlq_n_s16(cB_lo, 1)), rB_lo);
                int16x8_t by_hi = vaddq_s16(vaddq_s16(lB_hi, vshlq_n_s16(cB_hi, 1)), rB_hi);
                int16x8_t ty_lo = vaddq_s16(vaddq_s16(lT_lo, vshlq_n_s16(cT_lo, 1)), rT_lo);
                int16x8_t ty_hi = vaddq_s16(vaddq_s16(lT_hi, vshlq_n_s16(cT_hi, 1)), rT_hi);
                int16x8_t gy_lo = vsubq_s16(by_lo, ty_lo);
                int16x8_t gy_hi = vsubq_s16(by_hi, ty_hi);

                // L1 magnitude: |gx| + |gy|
                int16x8_t mag_lo = vqaddq_s16(vabsq_s16(gx_lo), vabsq_s16(gy_lo));
                int16x8_t mag_hi = vqaddq_s16(vabsq_s16(gx_hi), vabsq_s16(gy_hi));

                // Saturating narrow to u8
                uint8x8_t out_lo = vqmovun_s16(mag_lo);
                uint8x8_t out_hi = vqmovun_s16(mag_hi);
                vst1q_u8(rowO + x, vcombine_u8(out_lo, out_hi));

                // advance previous blocks
                prevT = curT; prevC = curC; prevB = curB;
            }

            // tail (scalar)
            for (; x < W; ++x) {
                int tL = (rowT && x-1 >= 0) ? rowT[x-1] : 0;
                int tC = (rowT) ? rowT[x] : 0;
                int tR = (rowT && x+1 < W) ? rowT[x+1] : 0;
                int cL = (x-1 >= 0) ? rowC[x-1] : 0;
                int cR = (x+1 <  W) ? rowC[x+1] : 0;
                int bL = (rowB && x-1 >= 0) ? rowB[x-1] : 0;
                int bC = (rowB) ? rowB[x] : 0;
                int bR = (rowB && x+1 < W) ? rowB[x+1] : 0;

                int gx = (tR + (cR<<1) + bR) - (tL + (cL<<1) + bL);
                int gy = (bL + (bC<<1) + bR) - (tL + (tC<<1) + tR);
                int mag = (gx<0?-gx:gx) + (gy<0?-gy:gy);
                rowO[x] = clamp_u8_int(mag);
            }
        } else {
            // narrow rows (scalar)
            for (int x = 0; x < W; ++x) {
                int tL = (rowT && x-1 >= 0) ? rowT[x-1] : 0;
                int tC = (rowT) ? rowT[x] : 0;
                int tR = (rowT && x+1 < W) ? rowT[x+1] : 0;
                int cL = (x-1 >= 0) ? rowC[x-1] : 0;
                int cR = (x+1 <  W) ? rowC[x+1] : 0;
                int bL = (rowB && x-1 >= 0) ? rowB[x-1] : 0;
                int bC = (rowB) ? rowB[x] : 0;
                int bR = (rowB && x+1 < W) ? rowB[x+1] : 0;

                int gx = (tR + (cR<<1) + bR) - (tL + (cL<<1) + bL);
                int gy = (bL + (bC<<1) + bR) - (tL + (tC<<1) + tR);
                int mag = (gx<0?-gx:gx) + (gy<0?-gy:gy);
                rowO[x] = clamp_u8_int(mag);
            }
        }
    }

    return out;
}

int boundaryCheck(int index_x, int index_y, int width, int height) {
    if(index_x >= 0 && index_y >= 0 && index_x < width && index_y < height) return 1;
    else return 0;
}
// Algorithms End.


/*******************************************************************************/
//Read PPM image and return an image pointer
/**************************************************************************/
Image *ReadPNMImage(char *filename)
{
    char ch;
    int  maxval, Width, Height;
    int size, num,j;
    FILE *fp;
    Image *image;
    int num_comment_lines=0;
    
    
    image=(Image *)malloc(sizeof(Image));
    
    if((fp=fopen(filename,"rb")) == NULL){
        printf("Cannot open %s\n", filename);
        exit(0);
    }
    
    printf("Loading %s ...",filename);
    
    if (fscanf(fp, "P%c\n", &ch) != 1) {
        printf("File is not in ppm/pgm raw format; cannot read\n");
        exit(0);
    }
    if( ch != '6' && ch !='5') {
        printf("File is not in ppm/pgm raw format; cannot read\n");
        exit(0);
    }
    
    if(ch == '5')image->Type=GRAY; 
    else if(ch == '6')image->Type=COLOR;  
    ch = getc(fp);
    j=0;
    while (ch == '#')
    {
        image->comments[num_comment_lines][j]=ch;
        j++;
        do {
            ch = getc(fp);
            image->comments[num_comment_lines][j]=ch;
            j++;
        } while (ch != '\n'); 
        image->comments[num_comment_lines][j-1]='\0';
        j=0;
        num_comment_lines++;
        ch = getc(fp);            
    }
    
    if (!isdigit((int)ch)){
        printf("Cannot read header information from ppm file");
        exit(0);
    }
    
    ungetc(ch, fp);       
    
    fscanf(fp, "%d%d%d\n", &Width, &Height, &maxval);
    
    if(image->Type == GRAY)
        size          = Width * Height;
    else  if(image->Type == COLOR)
        size          = Width * Height *3;
    image->data   = (unsigned char *) malloc(size);
    image->Width  = Width;
    image->Height = Height;
    image->num_comment_lines= num_comment_lines;
    
    if (!image->data){
        printf("cannot allocate memory for new image");
        exit(0);
    }
    
    num = fread((void *) image->data, 1, (size_t) size, fp);
    if (num != size){
        printf("cannot read image data from file");
        // exit(0);
    }
    
    fclose(fp);
    
    if(image->Type == GRAY)printf("..Image Type PGM\n");
    else printf("..Image Type PPM Color\n");
    return(image);
}

char *Extract_Filename_Stem(char *file_in) {
    const char *fname = file_in;
    const char *slash1 = strrchr(file_in, '/');
    const char *slash2 = strrchr(file_in, '\\'); 
    if (slash1 || slash2) {
        const char *slash = (slash1 > slash2 ? slash1 : slash2);
        fname = slash + 1;
    }

    char *stem = strdup(fname);
    if (!stem) return NULL;

    char *dot = strrchr(stem, '.');
    if (dot) *dot = '\0';

    return stem;
}

void SavePNMImage(Image *temp_image, char *filename, char *pgm_directory_name, char *png_directory_name)
{
    int num, j;
    int size;
    FILE *fp;

    ensure_output_dir(pgm_directory_name);
    ensure_output_dir(png_directory_name);
    char pgm_path[512];
    snprintf(pgm_path, sizeof(pgm_path), "%s/%s", pgm_directory_name, filename);

    const char *dot = strrchr(filename, '.');
    char png_filename[256];
    if (dot && (strcasecmp(dot, ".pgm") == 0 || strcasecmp(dot, ".ppm") == 0)) {
        size_t base_len = (size_t)(dot - filename);
        if (base_len >= sizeof(png_filename) - 5) base_len = sizeof(png_filename) - 5;
        memcpy(png_filename, filename, base_len);
        png_filename[base_len] = '\0';
        strncat(png_filename, ".png", sizeof(png_filename) - strlen(png_filename) - 1);
    } else {
        snprintf(png_filename, sizeof(png_filename), "%s.png", filename);
    }

    char png_path[512];
    snprintf(png_path, sizeof(png_path), "%s/%s", png_directory_name, png_filename);

    printf("Saving Image %s\n", pgm_path);
    fp = fopen(pgm_path, "wb");
    if (!fp) {
        perror("cannot open file for writing");
        exit(1);
    }

    if (temp_image->Type == GRAY) { 
        fprintf(fp, "P5\n");
        size = temp_image->Width * temp_image->Height;
    } else if (temp_image->Type == COLOR) { 
        fprintf(fp, "P6\n");
        size = temp_image->Width * temp_image->Height * 3;
    } else {
        fprintf(stderr, "Unsupported image type\n");
        fclose(fp);
        exit(1);
    }

    for (j = 0; j < temp_image->num_comment_lines; j++)
        fprintf(fp, "%s\n", temp_image->comments[j]);

    fprintf(fp, "%d %d\n%d\n", temp_image->Width, temp_image->Height, 255);

    num = (int)fwrite((void *)temp_image->data, 1, (size_t)size, fp);
    if (num != size) {
        fprintf(stderr, "cannot write image data to file\n");
        fclose(fp);
        exit(1);
    }
    fclose(fp);

    char cmd[1200];
    snprintf(cmd, sizeof(cmd), "magick \"%s\" \"%s\"", pgm_path, png_path);
    int rc = system(cmd);
    if (rc != 0) {
        fprintf(stderr, "ImageMagick conversion failed (code %d) for %s -> %s\n", rc, pgm_path, png_path);
    } else {
        printf("Saved PNG %s\n", png_path);
    }
}

static inline void neon_memset_u8(void *dst, uint8_t val, size_t n) {
    if (n == 0) return;

    uint8x16_t v = vdupq_n_u8(val);

    uint8_t *p = (uint8_t*)dst;
    while (n >= 64) {
        vst1q_u8(p +  0, v);
        vst1q_u8(p + 16, v);
        vst1q_u8(p + 32, v);
        vst1q_u8(p + 48, v);
        p += 64;
        n -= 64;
    }
    while (n >= 16) {
        vst1q_u8(p, v);
        p += 16;
        n -= 16;
    }
    if (n) {
        memset(p, val, n);
    }
}

/*************************************************************************/
/*Create a New Image with same dimensions as input image*/
/*************************************************************************/

Image *CreateNewImage(Image * image, char *comment)
{
    Image *outimage = (Image *)malloc(sizeof(Image));
    if (!outimage) {
        fprintf(stderr, "cannot allocate Image struct\n");
        exit(1);
    }

    outimage->Type = image->Type;
    outimage->num_comment_lines = image->num_comment_lines;

    int size = 0;
    if (outimage->Type == GRAY)  size = image->Width * image->Height;
    if (outimage->Type == COLOR) size = image->Width * image->Height * 3;

    outimage->Width  = image->Width;
    outimage->Height = image->Height;

    // Copy prior comments (strings are short; memcpy/strcpy is already optimal)
    for (int j = 0; j < outimage->num_comment_lines; j++) {
        strcpy(outimage->comments[j], image->comments[j]);
    }
    // Append new comment
    strcpy(outimage->comments[outimage->num_comment_lines], comment);
    outimage->num_comment_lines++;

    // 64-byte aligned allocation helps cache and wide stores
    void *buf = NULL;
    if (posix_memalign(&buf, 64, (size_t)size) != 0 || !buf) {
        fprintf(stderr, "cannot allocate memory for new image data\n");
        free(outimage);
        exit(1);
    }
    outimage->data = (unsigned char*)buf;

    // Fast zero-init with NEON (most predictable “SIMD win” in this function)
    neon_memset_u8(outimage->data, 0u, (size_t)size);

    return outimage;
}