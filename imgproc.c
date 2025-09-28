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

Image *Laplacian(Image *);
Image *Sobel(Image *);
Image *Gamma(Image *, float ratio);
Image *Global_histogram(Image *);
Image *Local_histogram(Image *);
Image *ReadPNMImage(char *);
Image *CreateNewImage(Image *, char *comment);
int boundaryCheck(int index_x, int index_y, int width, int height);
int TestReadImage(char *, char *);
void SavePNMImage(Image *, char *, char*, char*);
char *Extract_Filename_Stem(char*);
void process_folder_to_pgm_then_run(char *input_dir, char *pgm_output_dir);

int main(int argc, char **argv)
{
    // Please adjust the input filename and path to suit your needs:
    // char* file_in = (char*)"lena.pgm";
    // char* file_out = (char*)"";
    char *input_dir = (argc > 1) ? argv[1] : "images";
    char *pgm_output_dir = (argc > 2) ? argv[2] : "grayscale_inputs_pgm";
    process_folder_to_pgm_then_run(input_dir, pgm_output_dir);

    // TestReadImage(file_in, file_out);
    return(0);
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
        // skip . and ..
        if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0) continue;

        // Build full input path: <input_dir>/<name>
        char in_path[1024];
        snprintf(in_path, sizeof(in_path), "%s/%s", input_dir, ent->d_name);

        // Convert to PGM
        char pgm_path[1024];
        if (convert_to_pgm(in_path, pgm_output_dir, pgm_path, sizeof(pgm_path)) == 0) {
            // Now run your pipeline on the generated PGM
            // file_out is unused in your current TestReadImage; pass empty string
            TestReadImage(pgm_path, "");
        }
    }

    closedir(dir);
}

int TestReadImage(char *file_in, char *file_out)
{
    Image *image;
    Image *laplacian;
    // Image *sobel, *histogram_global, *histogram_local;
    // Image *gama_01, *gama_04, *gama_07, *gama_1;
    char *stem = Extract_Filename_Stem(file_in);
    image = ReadPNMImage(file_in);
    laplacian = Laplacian(image);
    // sobel = Sobel(image);
    // gama_01 = Gamma(image, 0.1);
    // gama_04 = Gamma(image, 0.4);
    // gama_07 = Gamma(image, 0.7);
    // gama_1 = Gamma(image, 1);
    // histogram_global = Global_histogram(image);
    // histogram_local = Local_histogram(image);

    // Please adjust the output filenames and paths to suit your needs:
    char filename[300];
    snprintf(filename, sizeof(filename), "%s.pgm", stem);

    SavePNMImage(laplacian, filename, "laplacian_pgm", "laplacian_png");
    // SavePNMImage(sobel, (char*)"sobel.pgm");
    // SavePNMImage(gama_01, (char*)"gama_01.pgm");
    // SavePNMImage(gama_04, (char*)"gama_04.pgm");
    // SavePNMImage(gama_07, (char*)"gama_07.pgm");
    // SavePNMImage(gama_1, (char*)"gama_1.pgm");
    // SavePNMImage(histogram_global, (char*)"histogram_global.pgm");
    // SavePNMImage(histogram_local, (char*)"histogram_local.pgm");
    free(stem);
    return(0);
}

// Algorithms Code:
Image *Laplacian(Image *image) {
    unsigned char *tempin, *tempout;
    int sum = 0;
    Image *outimage;
    outimage = CreateNewImage(image, (char*)"#testing function");
    tempin = image->data;
    tempout = outimage->data;
    
    for(int i = 0; i < image->Height; i++) {
        for(int j = 0; j < image->Width; j++) {
            sum = 0;
            for(int m = -1; m <= 1; m += 2) {
                for(int n = -1; n <= 1; n += 2) {
                    // use boundary check:
                    sum += boundaryCheck(j + n, i + m, image->Width, image->Height) ? tempin[image->Width * (i + m) + (j + n)] : 0;
                }
            }
            int temp = tempin[image->Width * i + j] * 4 - sum;
            // handle excess values:
            if(temp > 255) temp = 255;
            if(temp < 0) temp = 0;
            tempout[image->Width * i + j] = temp;
        }
    }
    return (outimage);
}

Image *Sobel(Image *image) {
    unsigned char *tempin, *tempout;
    int index, square[9], temp1, temp2;
    Image *outimage;
    outimage = CreateNewImage(image, (char*)"#testing function");
    tempin = image->data;
    tempout = outimage->data;
    
    for(int i = 0; i < image->Height; i++) {
        for(int j = 0; j < image->Width; j++) {
            index = 0;
            // record the values in the 3x3 square:
            for(int m = -1; m <= 1; m++) {
                for(int n = -1; n <= 1; n++) {
                    // use boundary check:
                    square[index++] = boundaryCheck(j + n, i + m, image->Width, image->Height) ? tempin[image->Width * (i + m) + (j + n)] : 0;
                }
            }
            temp1 = abs(square[2] + 2*square[5] + square[8] - square[0] - 2*square[3] - square[6]);
            temp2 = abs(square[6] + 2*square[7] + square[8] - square[0] - 2*square[1] - square[2]);
            tempout[image->Width * i + j] = sqrt(pow(temp1, 2) + pow(temp2, 2));
        }
    }
    return (outimage);
}

Image *Gamma(Image *image, float ratio) {
    unsigned char *tempin, *tempout;
    float temp;
    float variance, average, sum = 0, N = image->Width * image->Height; // calculate variance
    
    Image *outimage;
    outimage = CreateNewImage(image, (char*)"#testing function");
    tempin = image->data;
    tempout = outimage->data;
    
    for(int i = 0; i < image->Height; i++) {
        for(int j = 0; j < image->Width; j++) {
            temp = ((float)tempin[image->Width * i + j] + 0.5) / 256; // normalized
            temp = pow(temp, ratio); // power the parameter
            temp = (int)(temp * 256 - 0.5); // denormalization
            tempout[outimage->Width * i + j] = (unsigned char)temp;
            sum += temp;
        }
    }
    
    // calculate & output the variance:
    average = sum / N;
    sum = 0;
    for(int i = 0; i < image->Height; i++) {
        for(int j = 0; j < image->Width; j++) {
            sum += pow(tempout[outimage->Width * i + j] - average, 2);
        }
    }
    variance = sum / N;
    printf("The variance of gamma value %.1f is: %.2f\n", ratio, variance);
    return (outimage);
}

Image *Global_histogram(Image *image) {
    unsigned char *tempin, *tempout, temp;
    int histogram_sum[256], histogram[256], currSum = 0; // used for statistics
    float constant = (float)255 / (float)(image->Width * image->Height);// (L-1)/(M*N)
    
    Image *outimage;
    outimage = CreateNewImage(image, (char*)"#testing function");
    tempin = image->data;
    tempout = outimage->data;
    
    // initialize the array:
    for(int i = 0; i < 256; i++) histogram[i] = 0;
    
    for(int i = 0; i < image->Height; i++) {
        for(int j = 0; j < image->Width; j++) {
            temp = tempin[image->Width * i + j];
            histogram[temp] += 1;
        }
    }
    for(int i = 0; i < 256; i++) {
        currSum += histogram[i];
        histogram_sum[i] = currSum;
    }
    
    // output the image:
    for(int i = 0; i < image->Height; i++) {
        for(int j = 0; j < image->Width; j++) {
            temp = tempin[image->Width * i + j];
            tempout[outimage->Width * i + j] = (int)(histogram_sum[temp] * constant);
        }
    }
    return (outimage);
}

Image *Local_histogram(Image *image) {
    unsigned char *tempin, *tempout, temp;
    int histogram_sum[256], histogram[256]; // used for statistics
    float constant = (float)255 / (float)9; // M*N changed to 9
    
    Image *outimage;
    outimage = CreateNewImage(image, (char*)"#testing function");
    tempin = image->data;
    tempout = outimage->data;
    
    // process all the pixels:
    for(int i = 0; i < image->Height; i++) {
        for(int j = 0; j < image->Width; j++) {
            // initialize the array:
            for(int k = 0; k < 256; k++) histogram[k] = 0;

            for(int x = -1; x <= 1; x++) {
                for(int y = -1; y <= 1; y++) {
                    // use boundary check:
                    temp = boundaryCheck(j + y, i + x, image->Width, image->Height) ? tempin[image->Width * (i + x) + (j + y)] : 0;
                    histogram[temp] += 1;
                }
            }
            for(int k = 0, currSum = 0; k < 256; k++) {
                currSum += histogram[k];
                histogram_sum[k] = currSum;
            }
            
            // output the image:
            temp = tempin[image->Width * i + j];
            tempout[outimage->Width * i + j] = (int)(histogram_sum[temp] * constant);
        }
    }
    return (outimage);
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
    
    if(ch == '5')image->Type=GRAY;  // Gray (pgm)
    else if(ch == '6')image->Type=COLOR;  //Color (ppm)
    /* skip comments */
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
        } while (ch != '\n');     /* read to the end of the line */
        image->comments[num_comment_lines][j-1]='\0';
        j=0;
        num_comment_lines++;
        ch = getc(fp);            /* thanks, Elliot */
    }
    
    if (!isdigit((int)ch)){
        printf("Cannot read header information from ppm file");
        exit(0);
    }
    
    ungetc(ch, fp);               /* put that digit back */
    
    /* read the width, height, and maximum value for a pixel */
    fscanf(fp, "%d%d%d\n", &Width, &Height, &maxval);
    
    /*
     if (maxval != 255){
     printf("image is not true-color (24 bit); read failed");
     exit(0);
     }
     */
    
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
    //printf("Complete reading of %d bytes \n", num);
    if (num != size){
        printf("cannot read image data from file");
        exit(0);
    }
    
    //for(j=0;j<image->num_comment_lines;j++){
    //      printf("%s\n",image->comments[j]);
    //      }
    
    fclose(fp);
    
    /*-----  Debug  ------*/
    
    if(image->Type == GRAY)printf("..Image Type PGM\n");
    else printf("..Image Type PPM Color\n");
    /*
     printf("Width %d\n", Width);
     printf("Height %d\n",Height);
     printf("Size of image %d bytes\n",size);
     printf("maxvalue %d\n", maxval);
     */
    return(image);
}

char *Extract_Filename_Stem(char *file_in) {
    const char *fname = file_in;
    const char *slash1 = strrchr(file_in, '/');
    const char *slash2 = strrchr(file_in, '\\'); // Windows support
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

    // Ensure output directories exist
    ensure_output_dir(pgm_directory_name);
    ensure_output_dir(png_directory_name);

    // Build PGM full path: <pgm_directory_name>/<filename>
    char pgm_path[512];
    snprintf(pgm_path, sizeof(pgm_path), "%s/%s", pgm_directory_name, filename);

    // Derive PNG filename from `filename`
    // If filename ends with .pgm or .ppm -> replace with .png, else append .png
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

    // Build PNG full path: <png_directory_name>/<png_filename>
    char png_path[512];
    snprintf(png_path, sizeof(png_path), "%s/%s", png_directory_name, png_filename);

    // ---- Write PNM (PGM/PPM) ----
    printf("Saving Image %s\n", pgm_path);
    fp = fopen(pgm_path, "wb");
    if (!fp) {
        perror("cannot open file for writing");
        exit(1);
    }

    if (temp_image->Type == GRAY) {  // Gray (pgm)
        fprintf(fp, "P5\n");
        size = temp_image->Width * temp_image->Height;
    } else if (temp_image->Type == COLOR) {  // Color (ppm)
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

    // ---- Auto-convert PGM/PPM -> PNG via ImageMagick so you don't have to run it manually ----
    // Uses `magick` CLI; ensure ImageMagick is installed. Quotes handle spaces in paths.
    char cmd[1200];
    snprintf(cmd, sizeof(cmd), "magick \"%s\" \"%s\"", pgm_path, png_path);
    int rc = system(cmd);
    if (rc != 0) {
        fprintf(stderr, "ImageMagick conversion failed (code %d) for %s -> %s\n", rc, pgm_path, png_path);
    } else {
        printf("Saved PNG %s\n", png_path);
    }
}

/*************************************************************************/
/*Create a New Image with same dimensions as input image*/
/*************************************************************************/

Image *CreateNewImage(Image * image, char *comment)
{
    Image *outimage;
    int size,j;
    
    outimage=(Image *)malloc(sizeof(Image));
    
    outimage->Type = image->Type;
    outimage->num_comment_lines = image->num_comment_lines;

    if(outimage->Type == GRAY)   size = image->Width * image->Height;
    if(outimage->Type == COLOR) size  = image->Width * image->Height * 3;
    outimage->Width = image->Width;
    outimage->Height = image->Height;
    
    /*--------------------------------------------------------*/
    /* Copy Comments for Original Image      */
    for(j=0;j<outimage->num_comment_lines;j++)
        strcpy(outimage->comments[j],image->comments[j]);
    
    /*----------- Add New Comment  ---------------------------*/
    strcpy(outimage->comments[outimage->num_comment_lines],comment);
    outimage->num_comment_lines++;
    
    
    outimage->data = (unsigned char *) malloc(size);
    if (!outimage->data){
        printf("cannot allocate memory for new image");
        exit(0);
    }
    return(outimage);
}