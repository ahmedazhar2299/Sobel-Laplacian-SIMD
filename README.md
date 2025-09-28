# Sobel-Laplacian-SIMD

This project provides a simple image processing tool implementing **Sobel** and **Laplacian** edge detection filters using SIMD optimizations. It demonstrates how to preprocess images, compile C-based image processing code, and convert results into viewable formats.

---

## Installation

```bash
brew install imagemagick

# Compile the image processing source code
clang -O3 -std=c11 imgproc.c -o imgproc

# Run the program
./imgproc
```

## References

- [Wilson-ZheLin / Introduction-to-Digital-Image-Processing](https://github.com/Wilson-ZheLin/Introduction-to-Digital-Image-Processing/tree/main/4.%20Edge%20Detection%20and%20Grayscale%20Transformation)
- [ImageMagick Documentation](https://imagemagick.org/)
