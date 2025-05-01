#include <vector>      // For std::vector
#include <cmath>       // For fabs, log10
#include <iostream>    // For std::cout, std::cerr
#include <cstdint>     // For uint8_t

// stb_image and stb_image_write (define implementation in one source file)
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// Define macros for min/max
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

// Small epsilon to prevent division by zero
const float EPSILON = 1e-6f;

// Struct for RGB image (Section 1: Full-color output)
struct Image {
    std::vector<std::vector<float>> r; // Red channel
    std::vector<std::vector<float>> g; // Green channel
    std::vector<std::vector<float>> b; // Blue channel
    int width;                         // Image width
    int height;                        // Image height
};

// Struct for Bayer pattern (Section 1: Mosaic input)
struct BayerPattern {
    std::vector<std::vector<float>> data; // Single color values (R, G, or B)
    int width;                           // Pattern width
    int height;                          // Pattern height
};

// Create Image and initialize to zero
Image* create_image(int height, int width) {
    Image* img = new Image;
    img->width = width;
    img->height = height;
    img->r = std::vector<std::vector<float>>(height, std::vector<float>(width, 0.0f));
    img->g = std::vector<std::vector<float>>(height, std::vector<float>(width, 0.0f));
    img->b = std::vector<std::vector<float>>(height, std::vector<float>(width, 0.0f));
    return img;
}

// Free Image
void free_image(Image* img) {
    if (img) {
        delete img; // std::vector cleans up automatically
    }
}

// Create BayerPattern and initialize to zero
BayerPattern* create_bayer_pattern(int height, int width) {
    BayerPattern* bayer = new BayerPattern;
    bayer->width = width;
    bayer->height = height;
    bayer->data = std::vector<std::vector<float>>(height, std::vector<float>(width, 0.0f));
    return bayer;
}

// Free BayerPattern
void free_bayer_pattern(BayerPattern* bayer) {
    if (bayer) {
        delete bayer; // std::vector cleans up automatically
    }
}

// Get pixel value with boundary check (Section 2.2.2.1)
float get_value(int i, int j, char color, const Image* img) {
    if (i < 0 || i >= img->height || j < 0 || j >= img->width) return 0.0f;
    switch (color) {
        case 'r': return img->r[i][j];
        case 'g': return img->g[i][j];
        case 'b': return img->b[i][j];
        default: return 0.0f;
    }
}

// Compute simple gradient (Section 2.2.2.1)
float compute_gradient(int i, int j, int k, int m, char color, const Image* img) {
    return fabs(get_value(i, j, color, img) - get_value(k, m, color, img));
}

// Compute interpolation weight (Section 2.2.2.1)
float compute_weight(float gradient) {
    return 1.0f / (1.0f + gradient);
}

// Interpolate green pixels (Section 2.2.2.1)
void interpolate_green(Image* img, const BayerPattern* bayer) {
    for (int i = 0; i < img->height; i++) {
        for (int j = 0; j < img->width; j++) {
            if ((i % 2 == 0 && j % 2 == 1) || (i % 2 == 1 && j % 2 == 0)) {
                img->g[i][j] = bayer->data[i][j];
            } else {
                float g2 = get_value(i - 1, j, 'g', img); // Top
                float g4 = get_value(i, j - 1, 'g', img); // Left
                float g6 = get_value(i, j + 1, 'g', img); // Right
                float g8 = get_value(i + 1, j, 'g', img); // Bottom
                float e2 = compute_weight(compute_gradient(i, j, i - 1, j, 'g', img));
                float e4 = compute_weight(compute_gradient(i, j, i, j - 1, 'g', img));
                float e6 = compute_weight(compute_gradient(i, j, i, j + 1, 'g', img));
                float e8 = compute_weight(compute_gradient(i, j, i + 1, j, 'g', img));
                float sum_weights = e2 + e4 + e6 + e8;
                if (sum_weights > 0) {
                    img->g[i][j] = (e2 * g2 + e4 * g4 + e6 * g6 + e8 * g8) / sum_weights;
                } else {
                    // Fallback: Average available neighbors, avoid zero
                    int count = 0;
                    float sum = 0.0f;
                    if (g2 > 0) { sum += g2; count++; }
                    if (g4 > 0) { sum += g4; count++; }
                    if (g6 > 0) { sum += g6; count++; }
                    if (g8 > 0) { sum += g8; count++; }
                    img->g[i][j] = (count > 0) ? (sum / count) : EPSILON;
                }
                // Clamp to prevent invalid values
                img->g[i][j] = std::min(std::max(img->g[i][j], EPSILON), 1.0f);
            }
        }
    }
}

// Interpolate red and blue pixels (Section 2.2.2.2)
void interpolate_red_blue(Image* img, const BayerPattern* bayer) {
    for (int i = 0; i < img->height; i++) {
        for (int j = 0; j < img->width; j++) {
            if (i % 2 == 0 && j % 2 == 0) {
                img->r[i][j] = bayer->data[i][j];
                img->b[i][j] = 0.0f;
            } else if (i % 2 == 1 && j % 2 == 1) {
                img->b[i][j] = bayer->data[i][j];
                img->r[i][j] = 0.0f;
            } else {
                img->r[i][j] = 0.0f;
                img->b[i][j] = 0.0f;
            }
            if (img->r[i][j] == 0.0f) {
                float sum_r_g = 0.0f, sum_weights = 0.0f;
                for (int di = -1; di <= 1; di += 2) {
                    for (int dj = -1; dj <= 1; dj += 2) {
                        float r = get_value(i + di, j + dj, 'r', img);
                        float g = get_value(i + di, j + dj, 'g', img);
                        if (g > EPSILON) { // Avoid division by zero
                            float weight = compute_weight(compute_gradient(i, j, i + di, j + dj, 'g', img));
                            sum_r_g += weight * (r / g);
                            sum_weights += weight;
                        }
                    }
                }
                if (sum_weights > 0) {
                    img->r[i][j] = img->g[i][j] * (sum_r_g / sum_weights);
                } else {
                    img->r[i][j] = img->g[i][j]; // Fallback: Use green value
                }
                img->r[i][j] = std::min(std::max(img->r[i][j], 0.0f), 1.0f);
            }
            if (img->b[i][j] == 0.0f) {
                float sum_b_g = 0.0f, sum_weights = 0.0f;
                for (int di = -1; di <= 1; di += 2) {
                    for (int dj = -1; dj <= 1; dj += 2) {
                        float b = get_value(i + di, j + dj, 'b', img);
                        float g = get_value(i + di, j + dj, 'g', img);
                        if (g > EPSILON) { // Avoid division by zero
                            float weight = compute_weight(compute_gradient(i, j, i + di, j + dj, 'g', img));
                            sum_b_g += weight * (b / g);
                            sum_weights += weight;
                        }
                    }
                }
                if (sum_weights > 0) {
                    img->b[i][j] = img->g[i][j] * (sum_b_g / sum_weights);
                } else {
                    img->b[i][j] = img->g[i][j]; // Fallback: Use green value
                }
                img->b[i][j] = std::min(std::max(img->b[i][j], 0.0f), 1.0f);
            }
        }
    }
}

// Correct colors with fixed iterations (Section 2.2.2.3)
void correct_colors(Image* img) {
    for (int iter = 0; iter < 3; iter++) {
        Image* temp = create_image(img->height, img->width);
        if (!temp) return;
        temp->r = img->r;
        temp->g = img->g;
        temp->b = img->b;
        for (int i = 0; i < img->height; i++) {
            for (int j = 0; j < img->width; j++) {
                float g_r = 0.0f, g_b = 0.0f, sum_weights = 0.0f;
                for (int di = -1; di <= 1; di += 2) {
                    for (int dj = -1; dj <= 1; dj += 2) {
                        float r = get_value(i + di, j + dj, 'r', temp);
                        float g = get_value(i + di, j + dj, 'g', temp);
                        float b = get_value(i + di, j + dj, 'b', temp);
                        if (g > EPSILON && r > EPSILON && b > EPSILON) { // Avoid division by zero
                            float weight = compute_weight(compute_gradient(i, j, i + di, j + dj, 'g', temp));
                            g_r += weight * (g / r);
                            g_b += weight * (g / b);
                            sum_weights += weight;
                        }
                    }
                }
                if (sum_weights > 0) {
                    float r5 = get_value(i, j, 'r', temp);
                    float b5 = get_value(i, j, 'b', temp);
                    float g_r_new = r5 * (g_r / sum_weights);
                    float g_b_new = b5 * (g_b / sum_weights);
                    img->g[i][j] = (g_r_new + g_b_new) / 2.0f;
                    img->g[i][j] = std::min(std::max(img->g[i][j], EPSILON), 1.0f);
                }
                float sum_r_g = 0.0f, sum_b_g = 0.0f;
                sum_weights = 0.0f;
                for (int di = -1; di <= 1; di += 2) {
                    for (int dj = -1; dj <= 1; dj += 2) {
                        float r = get_value(i + di, j + dj, 'r', temp);
                        float g = get_value(i + di, j + dj, 'g', temp);
                        float b = get_value(i + di, j + dj, 'b', temp);
                        if (g > EPSILON) { // Avoid division by zero
                            float weight = compute_weight(compute_gradient(i, j, i + di, j + dj, 'g', temp));
                            sum_r_g += weight * (r / g);
                            sum_b_g += weight * (b / g);
                            sum_weights += weight;
                        }
                    }
                }
                if (sum_weights > 0) {
                    img->r[i][j] = img->g[i][j] * (sum_r_g / sum_weights);
                    img->b[i][j] = img->g[i][j] * (sum_b_g / sum_weights);
                } else {
                    // Fallback: Retain current values
                    img->r[i][j] = img->r[i][j] > 0 ? img->r[i][j] : img->g[i][j];
                    img->b[i][j] = img->b[i][j] > 0 ? img->b[i][j] : img->g[i][j];
                }
                img->r[i][j] = std::min(std::max(img->r[i][j], 0.0f), 1.0f);
                img->b[i][j] = std::min(std::max(img->b[i][j], 0.0f), 1.0f);
            }
        }
        free_image(temp);
    }
}

// Main algorithm function (Section 2.2.2)
void kimmel_algorithm(Image* img, const BayerPattern* bayer) {
    interpolate_green(img, bayer);
    interpolate_red_blue(img, bayer);
    correct_colors(img);
}

// Convert Image to uint8_t array for saving
void image_to_buffer(const Image* img, uint8_t* buffer) {
    for (int i = 0; i < img->height; i++) {
        for (int j = 0; j < img->width; j++) {
            int idx = (i * img->width + j) * 3;
            // Convert float [0,1] to uint8 [0,255], RGB format
            buffer[idx]     = (uint8_t)(std::min(std::max(img->r[i][j], 0.0f), 1.0f) * 255.0f); // R
            buffer[idx + 1] = (uint8_t)(std::min(std::max(img->g[i][j], 0.0f), 1.0f) * 255.0f); // G
            buffer[idx + 2] = (uint8_t)(std::min(std::max(img->b[i][j], 0.0f), 1.0f) * 255.0f); // B
        }
    }
}

// Create Bayer pattern from color or grayscale image
void create_bayer_pattern(const uint8_t* data, int width, int height, int channels, BayerPattern* bayer) {
    bool is_grayscale = (channels == 1);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            int idx = (i * width + j) * channels;
            if (is_grayscale) {
                float intensity = data[idx] / 255.0f; // Normalize to [0,1]
                // GRBG pattern with adjusted color variation
                float r_scale = 0.9f; // Balanced red
                float g_scale = 1.1f; // Slightly emphasize green
                float b_scale = 0.95f; // Balanced blue
                if (i % 2 == 0 && j % 2 == 0) bayer->data[i][j] = intensity * r_scale; // Red
                else if (i % 2 == 0 && j % 2 == 1) bayer->data[i][j] = intensity * g_scale; // Green
                else if (i % 2 == 1 && j % 2 == 0) bayer->data[i][j] = intensity * g_scale; // Green
                else bayer->data[i][j] = intensity * b_scale; // Blue
            } else {
                float r = data[idx] / 255.0f;     // R (stb_image uses RGB)
                float g = data[idx + 1] / 255.0f; // G
                float b = data[idx + 2] / 255.0f; // B
                // GRBG pattern
                if (i % 2 == 0 && j % 2 == 0) bayer->data[i][j] = r; // Red
                else if (i % 2 == 0 && j % 2 == 1) bayer->data[i][j] = g; // Green
                else if (i % 2 == 1 && j % 2 == 0) bayer->data[i][j] = g; // Green
                else bayer->data[i][j] = b; // Blue
            }
            // Ensure non-zero values
            bayer->data[i][j] = std::max(bayer->data[i][j], EPSILON);
        }
    }
}

// Compute PSNR (color or grayscale)
double compute_psnr(const uint8_t* input, int input_channels, const Image* output, int width, int height) {
    double mse = 0.0;
    int count = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            int idx = (i * width + j) * input_channels;
            if (input_channels == 1) {
                // Grayscale: Compare input with output's grayscale (R=G=B)
                float input_val = input[idx] / 255.0f;
                float output_val = (output->r[i][j] + output->g[i][j] + output->b[i][j]) / 3.0f;
                if (std::isfinite(input_val) && std::isfinite(output_val)) {
                    float diff = input_val - output_val;
                    mse += diff * diff;
                    count++;
                }
            } else {
                // Color: Compare RGB
                float input_r = input[idx] / 255.0f;
                float input_g = input[idx + 1] / 255.0f;
                float input_b = input[idx + 2] / 255.0f;
                float diff_r = input_r - output->r[i][j];
                float diff_g = input_g - output->g[i][j];
                float diff_b = input_b - output->b[i][j];
                if (std::isfinite(diff_r) && std::isfinite(diff_g) && std::isfinite(diff_b)) {
                    mse += diff_r * diff_r + diff_g * diff_g + diff_b * diff_b;
                    count += 3;
                }
            }
        }
    }
    if (count == 0) return 0.0; // Avoid division by zero
    mse /= count;
    if (mse == 0) return 100.0; // Perfect match
    return 10.0 * log10((1.0 * 1.0) / mse); // Normalized to [0,1]
}

// Main function for testing with color or grayscale PNG (512x512)
int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_image.png>" << std::endl;
        return 1;
    }

    // Load image
    int width, height, channels;
    uint8_t* input_data = stbi_load(argv[1], &width, &height, &channels, 0);
    if (!input_data) {
        std::cerr << "Error: Could not load image " << argv[1] << std::endl;
        return 1;
    }
    if (width != 512 || height != 512) {
        std::cerr << "Error: Image must be 512x512" << std::endl;
        stbi_image_free(input_data);
        return 1;
    }
    if (channels != 1 && channels != 3) {
        std::cerr << "Error: Image must be grayscale (1 channel) or color (3 channels)" << std::endl;
        stbi_image_free(input_data);
        return 1;
    }

    // Create Bayer pattern
    BayerPattern* bayer = create_bayer_pattern(height, width);
    if (!bayer) {
        std::cerr << "Error: Failed to create Bayer pattern" << std::endl;
        stbi_image_free(input_data);
        return 1;
    }
    create_bayer_pattern(input_data, width, height, channels, bayer);

    // Save Bayer pattern for debugging
    uint8_t* bayer_buffer = new uint8_t[width * height];
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++)
            bayer_buffer[i * width + j] = (uint8_t)(bayer->data[i][j] * 255);
    stbi_write_png("bayer.png", width, height, 1, bayer_buffer, width);
    delete[] bayer_buffer;

    // Create output image
    Image* img = create_image(height, width);
    if (!img) {
        free_bayer_pattern(bayer);
        stbi_image_free(input_data);
        std::cerr << "Error: Failed to create image" << std::endl;
        return 1;
    }

    // Run demosaicing algorithm
    kimmel_algorithm(img, bayer);

    // Save output
    uint8_t* output_buffer = new uint8_t[width * height * 3];
    image_to_buffer(img, output_buffer);
    if (!stbi_write_png("output.png", width, height, 3, output_buffer, width * 3)) {
        std::cerr << "Error: Failed to save output.png" << std::endl;
        delete[] output_buffer;
        free_image(img);
        free_bayer_pattern(bayer);
        stbi_image_free(input_data);
        return 1;
    }
    delete[] output_buffer;

    // Compute PSNR
    double psnr = compute_psnr(input_data, channels, img, width, height);
    std::cout << "PSNR: " << psnr << " dB" << std::endl;

    // Print sample pixels
    std::cout << "Sample pixel (0,0): R=" << img->r[0][0]
              << ", G=" << img->g[0][0]
              << ", B=" << img->b[0][0] << std::endl;
    std::cout << "Sample pixel (1,1): R=" << img->r[1][1]
              << ", G=" << img->g[1][1]
              << ", B=" << img->b[1][1] << std::endl;

    // Clean up
    free_image(img);
    free_bayer_pattern(bayer);
    stbi_image_free(input_data);
    return 0;
}