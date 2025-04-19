#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdint>

struct Image {
    std::vector<std::vector<float>> r, g, b; // Color planes
    int width, height;
    Image(int w, int h) : width(w), height(h) {
        r.resize(h, std::vector<float>(w, 0.0f));
        g.resize(h, std::vector<float>(w, 0.0f));
        b.resize(h, std::vector<float>(w, 0.0f));
    }
};

struct BayerPattern {
    std::vector<std::vector<float>> data; // Bayer mosaic
    int width, height;
    BayerPattern(int w, int h) : width(w), height(h) {
        data.resize(h, std::vector<float>(w, 0.0f));
    }
};

// Compute gradients in a larger spatial window as per Section 3.2.2
float computeGradient(const Image& img, int i, int j, int k, int m, char color) {
    auto getValue = [&](int x, int y, char c) -> float {
        if (x < 0 || x >= img.height || y < 0 || y >= img.width) return 0.0f;
        if (c == 'G') return img.g[x][y];
        if (c == 'R') return img.r[x][y];
        return img.b[x][y];
    };
    
    float grad = 8.0f * std::abs(getValue(i, j, color) - getValue(k, m, color));
    grad += std::abs(getValue(i + 1, j, color) - getValue(k + 1, m, color));
    grad += std::abs(getValue(i - 1, j, color) - getValue(k - 1, m, color));
    grad += std::abs(getValue(i, j + 1, color) - getValue(k, m + 1, color));
    grad += std::abs(getValue(i, j - 1, color) - getValue(k, m - 1, color));
    return grad;
}

// Compute weights based on derivatives (Section 2.2.2.1)
float computeWeight(const Image& img, int i, int j, char color, float derivative) {
    float d_center = 0.0f;
    if (color == 'G' && i >= 0 && i < img.height && j >= 0 && j < img.width) {
        d_center = std::abs(img.g[i][j]);
    }
    return 1.0f / std::sqrt(1.0f + derivative * derivative + d_center * d_center);
}

// Interpolate green color (Section 2.2.2.1)
void interpolateGreen(Image& img, const BayerPattern& bayer) {
    for (int i = 1; i < img.height - 1; i += 2) {
        for (int j = 1; j < img.width - 1; j += 2) {
            if ((i + j) % 2 == 0) continue; // Skip green pixels
            float weights[4], values[4];
            int di[4] = {-1, 0, 1, 0}, dj[4] = {0, 1, 0, -1};
            float sum_weights = 0.0f, sum_values = 0.0f;
            
            for (int k = 0; k < 4; k++) {
                int ni = i + di[k], nj = j + dj[k];
                if (ni < 0 || ni >= img.height || nj < 0 || nj >= img.width) continue;
                float grad = computeGradient(img, i, j, ni, nj, 'G');
                weights[k] = computeWeight(img, ni, nj, 'G', grad);
                values[k] = img.g[ni][nj];
                sum_weights += weights[k];
                sum_values += weights[k] * values[k];
            }
            
            if (sum_weights > 0) {
                img.g[i][j] = sum_values / sum_weights;
            }
        }
    }
}

// Interpolate red and blue using green (Section 2.2.2.2)
void interpolateRedBlue(Image& img) {
    for (int i = 1; i < img.height - 1; i++) {
        for (int j = 1; j < img.width - 1; j++) {
            if (img.g[i][j] == 0) continue; // Skip if green not interpolated
            float weights[4], ratios[4];
            int di[4] = {-1, 0, 1, 0}, dj[4] = {0, 1, 0, -1};
            float sum_weights_r = 0.0f, sum_ratios_r = 0.0f;
            float sum_weights_b = 0.0f, sum_ratios_b = 0.0f;
            
            for (int k = 0; k < 4; k++) {
                int ni = i + di[k], nj = j + dj[k];
                if (ni < 0 || ni >= img.height || nj < 0 || nj >= img.width) continue;
                float grad_r = computeGradient(img, i, j, ni, nj, 'R');
                float grad_b = computeGradient(img, i, j, ni, nj, 'B');
                weights[k] = computeWeight(img, ni, nj, 'G', grad_r);
                
                if (img.g[ni][nj] > 0) {
                    if (img.r[ni][nj] > 0) {
                        ratios[k] = img.r[ni][nj] / img.g[ni][nj];
                        sum_ratios_r += weights[k] * ratios[k];
                        sum_weights_r += weights[k];
                    }
                    if (img.b[ni][nj] > 0) {
                        ratios[k] = img.b[ni][nj] / img.g[ni][nj];
                        sum_ratios_b += weights[k] * ratios[k];
                        sum_weights_b += weights[k];
                    }
                }
            }
            
            if (sum_weights_r > 0) img.r[i][j] = img.g[i][j] * (sum_ratios_r / sum_weights_r);
            if (sum_weights_b > 0) img.b[i][j] = img.g[i][j] * (sum_ratios_b / sum_weights_b);
        }
    }
}

// Correction stage with adaptive iterations (Section 2.2.2.3, 3.2.1)
void correctColors(Image& img, const std::vector<std::vector<float>>& problem_map) {
    Image temp(img.width, img.height);
    int max_iterations = 12;
    int min_iterations = 1;
    
    for (int i = 0; i < img.height; i++) {
        for (int j = 0; j < img.width; j++) {
            temp.r[i][j] = img.r[i][j];
            temp.g[i][j] = img.g[i][j];
            temp.b[i][j] = img.b[i][j];
        }
    }
    
    for (int i = 1; i < img.height - 1; i++) {
        for (int j = 1; j < img.width - 1; j++) {
            int iterations = min_iterations + static_cast<int>(problem_map[i][j] * (max_iterations - min_iterations));
            bool forward = true;
            
            for (int iter = 0; iter < iterations; iter++, forward = !forward) {
                int start_i = forward ? 1 : img.height - 2;
                int end_i = forward ? img.height - 1 : 1;
                int step_i = forward ? 1 : -1;
                
                for (int ii = start_i; forward ? ii < end_i : ii >= end_i; ii += step_i) {
                    for (int jj = (forward ? 1 : img.width - 2); 
                         forward ? jj < img.width - 1 : jj >= 1; 
                         jj += (forward ? 1 : -1)) {
                        float weights[4], g_r_ratios[4], g_b_ratios[4];
                        int di[4] = {-1, 0, 1, 0}, dj[4] = {0, 1, 0, -1};
                        float sum_weights = 0.0f, sum_g_r = 0.0f, sum_g_b = 0.0f;
                        
                        for (int k = 0; k < 4; k++) {
                            int ni = ii + di[k], nj = jj + dj[k];
                            if (ni < 0 || ni >= img.height || nj < 0 || nj >= img.width) continue;
                            float grad = computeGradient(img, ii, jj, ni, nj, 'G');
                            weights[k] = computeWeight(img, ni, nj, 'G', grad);
                            
                            if (img.r[ni][nj] > 0 && img.g[ni][nj] > 0) {
                                g_r_ratios[k] = img.g[ni][nj] / img.r[ni][nj];
                                sum_g_r += weights[k] * g_r_ratios[k];
                            }
                            if (img.b[ni][nj] > 0 && img.g[ni][nj] > 0) {
                                g_b_ratios[k] = img.g[ni][nj] / img.b[ni][nj];
                                sum_g_b += weights[k] * g_b_ratios[k];
                            }
                            sum_weights += weights[k];
                        }
                        
                        if (sum_weights > 0) {
                            float g_r = img.r[ii][jj] * (sum_g_r / sum_weights);
                            float g_b = img.b[ii][jj] * (sum_g_b / sum_weights);
                            img.g[ii][jj] = (g_r + g_b) / 2.0f;
                            
                            // Update red and blue
                            sum_weights = 0.0f;
                            float sum_r_ratios = 0.0f, sum_b_ratios = 0.0f;
                            for (int k = 0; k < 4; k++) {
                                int ni = ii + di[k], nj = jj + dj[k];
                                if (ni < 0 || ni >= img.height || nj < 0 || nj >= img.width) continue;
                                if (img.g[ni][nj] > 0) {
                                    if (img.r[ni][nj] > 0) {
                                        sum_r_ratios += weights[k] * (img.r[ni][nj] / img.g[ni][nj]);
                                    }
                                    if (img.b[ni][nj] > 0) {
                                        sum_b_ratios += weights[k] * (img.b[ni][nj] / img.g[ni][nj]);
                                    }
                                    sum_weights += weights[k];
                                }
                            }
                            
                            if (sum_weights > 0) {
                                img.r[ii][jj] = img.g[ii][jj] * (sum_r_ratios / sum_weights);
                                img.b[ii][jj] = img.g[ii][jj] * (sum_b_ratios / sum_weights);
                            }
                        }
                    }
                }
            }
        }
    }
}

// Modified Kimmel algorithm
void modifiedKimmelAlgorithm(Image& img, const BayerPattern& bayer, const std::vector<std::vector<float>>& problem_map) {
    // Step 1: Interpolate green
    interpolateGreen(img, bayer);
    
    // Step 2: Interpolate red and blue
    interpolateRedBlue(img);
    
    // Step 3: Correct colors with adaptive iterations
    correctColors(img, problem_map);
}

int main() {
    int width = 512, height = 512;
    BayerPattern bayer(width, height);
    Image img(width, height);
    std::vector<std::vector<float>> problem_map(height, std::vector<float>(width, 0.5f)); // Dummy problem map
    
    // Initialize bayer and img with data (not shown)
    modifiedKimmelAlgorithm(img, bayer, problem_map);
    
    return 0;
}