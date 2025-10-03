#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <random>

// Define M_PI if not defined (common issue on Windows)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Vector3D class for 3D calculations
struct Vector3D {
    double x, y, z;
    
    Vector3D() : x(0), y(0), z(0) {}
    Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}
    
    Vector3D operator+(const Vector3D& v) const { return Vector3D(x + v.x, y + v.y, z + v.z); }
    Vector3D operator-(const Vector3D& v) const { return Vector3D(x - v.x, y - v.y, z - v.z); }
    Vector3D operator*(double t) const { return Vector3D(x * t, y * t, z * t); }
    Vector3D operator/(double t) const { return Vector3D(x / t, y / t, z / t); }
    
    double dot(const Vector3D& v) const { return x * v.x + y * v.y + z * v.z; }
    double length() const { return sqrt(x * x + y * y + z * z); }
    Vector3D normalize() const { double len = length(); return len > 0 ? *this / len : Vector3D(); }
};

// Color structure for rendering
struct Color {
    double r, g, b;
    Color() : r(0), g(0), b(0) {}
    Color(double r, double g, double b) : r(r), g(g), b(b) {}
    
    Color operator+(const Color& c) const { return Color(r + c.r, g + c.g, b + c.b); }
    Color operator*(double t) const { return Color(r * t, g * t, b * t); }
    
    Color clamp() const {
        return Color(
            std::max(0.0, std::min(1.0, r)),
            std::max(0.0, std::min(1.0, g)),
            std::max(0.0, std::min(1.0, b))
        );
    }
};

// Ray structure for raytracing
struct Ray {
    Vector3D position;
    Vector3D direction;
    
    Ray(const Vector3D& pos, const Vector3D& dir) : position(pos), direction(dir) {}
};

// Black hole parameters
const double SCHWARZSCHILD_RADIUS = 6;  // 2GM/c^2 in natural units
const double BLACK_HOLE_MASS = 5.0;       // Mass in natural units
const double STEP_SIZE = 0.01;            // Integration step size
const int MAX_STEPS = 10000;              // Maximum ray tracing steps

// Camera and rendering parameters
const int IMAGE_WIDTH = 800;
const int IMAGE_HEIGHT = 600;
const double CAMERA_DISTANCE = 1000.0;      // Distance from black hole
const double FOV = 60.0;                  // Field of view in degrees

class BlackHoleVisualizer {
private:
    std::vector<std::vector<Color>> image;
    std::mt19937 rng;
    
    // Calculate Schwarzschild metric coefficients
    double getMetricCoefficient(double r) const {
        return 1.0 - SCHWARZSCHILD_RADIUS / r;
    }
    
    // Calculate gravitational acceleration (simplified)
    Vector3D getGravitationalAcceleration(const Vector3D& pos) const {
        double r = pos.length();
        if (r < SCHWARZSCHILD_RADIUS * 1.1) return Vector3D(); // Inside event horizon
        
        double factor = -BLACK_HOLE_MASS / (r * r * r);
        // Include first-order relativistic correction
        double correction = 1.0 + 1.5 * SCHWARZSCHILD_RADIUS / r;
        return pos * (factor * correction);
    }
    
    // Trace a light ray through curved spacetime
    bool traceRay(Ray& ray, double& distance) const {
        Vector3D pos = ray.position;
        Vector3D vel = ray.direction;
        distance = 0.0;
        
        for (int step = 0; step < MAX_STEPS; step++) {
            double r = pos.length();
            
            // Check if ray hit the black hole (event horizon)
            if (r < SCHWARZSCHILD_RADIUS * 1.05) {
                return false; // Ray absorbed by black hole
            }
            
            // Check if ray escaped to infinity
            if (r > 100.0) {
                ray.position = pos;
                ray.direction = vel.normalize();
                return true;
            }
            
            // Calculate gravitational acceleration
            Vector3D accel = getGravitationalAcceleration(pos);
            
            // Integrate using Verlet method for better stability
            Vector3D newPos = pos + vel * STEP_SIZE + accel * (0.5 * STEP_SIZE * STEP_SIZE);
            Vector3D newAccel = getGravitationalAcceleration(newPos);
            Vector3D newVel = vel + (accel + newAccel) * (0.5 * STEP_SIZE);
            
            // Update position and velocity
            pos = newPos;
            vel = newVel;
            distance += STEP_SIZE;
        }
        
        ray.position = pos;
        ray.direction = vel.normalize();
        return true;
    }
    
    // Generate background star field
    Color getBackgroundColor(const Vector3D& direction) const {
        // Create a procedural star field
        double u = atan2(direction.z, direction.x) / (2 * M_PI) + 0.5;
        double v = acos(direction.y) / M_PI;
        
        // Simple noise-based stars
        int x = (int)(u * 1000) % 1000;
        int y = (int)(v * 1000) % 1000;
        
        double noise = sin(x * 0.1) * cos(y * 0.1) + 
                      sin(x * 0.23) * cos(y * 0.17) * 0.5 +
                      sin(x * 0.41) * cos(y * 0.31) * 0.25;
        
        if (noise > 1.5) {
            double intensity = (noise - 1.5) * 2.0;
            return Color(intensity, intensity * 0.9, intensity * 0.8);
        }
        
        // Background color (deep space)
        return Color(0.02, 0.02, 0.05);
    }
    
    // Create accretion disk around black hole
    Color getAccretionDiskColor(const Vector3D& pos) const {
        double r = sqrt(pos.x * pos.x + pos.z * pos.z); // Distance from rotation axis
        double height = abs(pos.y);
        
        // Accretion disk extends from 3 to 20 Schwarzschild radii
        double innerRadius = 3.0 * SCHWARZSCHILD_RADIUS;
        double outerRadius = 20.0 * SCHWARZSCHILD_RADIUS;
        double diskHeight = 0.5 * SCHWARZSCHILD_RADIUS;
        
        if (r > innerRadius && r < outerRadius && height < diskHeight) {
            // Temperature decreases with radius (T ∝ r^-3/4)
            double temp = pow(innerRadius / r, 0.75);
            
            // Color based on temperature (blackbody radiation)
            if (temp > 0.8) {
                return Color(1.0, 0.9, 0.7) * (temp * temp); // Hot white/blue
            } else if (temp > 0.4) {
                return Color(1.0, 0.7, 0.3) * (temp * 1.5); // Orange
            } else {
                return Color(0.8, 0.3, 0.1) * (temp * 2.0); // Red
            }
        }
        
        return Color(); // No disk contribution
    }
    
public:
    BlackHoleVisualizer() : rng(std::random_device{}()) {
        image.resize(IMAGE_HEIGHT, std::vector<Color>(IMAGE_WIDTH));
    }
    
    void render() {
        std::cout << "Rendering black hole visualization..." << std::endl;
        
        // Camera setup
        Vector3D cameraPos(0, 0, CAMERA_DISTANCE);
        Vector3D cameraDir(0, 0, -1);
        Vector3D cameraUp(0, 1, 0);
        Vector3D cameraRight = cameraDir.normalize().dot(cameraUp) < 0.99 ? 
                               Vector3D(1, 0, 0) : Vector3D(0, 1, 0);
        
        double aspectRatio = (double)IMAGE_WIDTH / IMAGE_HEIGHT;
        double fovRad = FOV * M_PI / 180.0;
        double viewportHeight = 2.0 * tan(fovRad / 2.0);
        double viewportWidth = aspectRatio * viewportHeight;
        
        for (int y = 0; y < IMAGE_HEIGHT; y++) {
            if (y % 50 == 0) {
                std::cout << "Progress: " << (100 * y / IMAGE_HEIGHT) << "%" << std::endl;
            }
            
            for (int x = 0; x < IMAGE_WIDTH; x++) {
                // Calculate ray direction
                double u = ((double)x / IMAGE_WIDTH - 0.5) * viewportWidth;
                double v = ((double)y / IMAGE_HEIGHT - 0.5) * viewportHeight;
                
                Vector3D rayDir = (cameraDir + cameraRight * u + cameraUp * v).normalize();
                Ray ray(cameraPos, rayDir);
                
                Color pixelColor;
                double distance;
                
                if (traceRay(ray, distance)) {
                    // Ray escaped - sample background or accretion disk
                    Color diskColor = getAccretionDiskColor(ray.position);
                    Color bgColor = getBackgroundColor(ray.direction);
                    
                    // Gravitational lensing effect - multiple images
                    if (distance > 5.0 && ray.position.length() < 50.0) {
                        // Add gravitational lensing brightening
                        double lensingFactor = 1.0 + 2.0 * SCHWARZSCHILD_RADIUS / ray.position.length();
                        bgColor = bgColor * lensingFactor;
                    }
                    
                    pixelColor = diskColor + bgColor;
                } else {
                    // Ray hit black hole - pure black with slight edge glow
                    double edgeGlow = exp(-distance * 10.0);
                    pixelColor = Color(edgeGlow * 0.1, edgeGlow * 0.05, 0.0);
                }
                
                // Apply relativistic Doppler shift effect (simplified)
                double dopplerShift = 1.0 - 0.1 * (u * u + v * v);
                if (dopplerShift < 0.7) dopplerShift = 0.7;
                pixelColor = Color(
                    pixelColor.r * dopplerShift,
                    pixelColor.g,
                    pixelColor.b / dopplerShift
                );
                
                image[y][x] = pixelColor.clamp();
            }
        }
        
        std::cout << "Rendering complete!" << std::endl;
    }
    
    void saveToPPM(const std::string& filename) const {
        std::ofstream file(filename);
        file << "P3\n" << IMAGE_WIDTH << " " << IMAGE_HEIGHT << "\n255\n";
        
        for (int y = 0; y < IMAGE_HEIGHT; y++) {
            for (int x = 0; x < IMAGE_WIDTH; x++) {
                Color c = image[y][x];
                int r = (int)(c.r * 255);
                int g = (int)(c.g * 255);
                int b = (int)(c.b * 255);
                file << r << " " << g << " " << b << "\n";
            }
        }
        
        file.close();
        std::cout << "Image saved to " << filename << std::endl;
    }
    
    // Render spacetime curvature visualization
    void renderSpacetimeCurvature(const std::string& filename) const {
        std::ofstream file(filename + "_curvature.txt");
        file << "# Spacetime curvature data\n";
        file << "# Format: x y z curvature\n";
        
        // Sample spacetime curvature on a grid
        for (double x = -10; x <= 10; x += 0.5) {
            for (double y = -10; y <= 10; y += 0.5) {
                for (double z = -10; z <= 10; z += 0.5) {
                    Vector3D pos(x, y, z);
                    double r = pos.length();
                    
                    if (r > SCHWARZSCHILD_RADIUS && r < 20.0) {
                        // Simplified curvature measure
                        double curvature = SCHWARZSCHILD_RADIUS / (r * r);
                        file << x << " " << y << " " << z << " " << curvature << "\n";
                    }
                }
            }
        }
        
        file.close();
        std::cout << "Spacetime curvature data saved to " << filename << "_curvature.txt" << std::endl;
    }
};

int main() {
    std::cout << "Black Hole Spacetime Visualization\n";
    std::cout << "==================================\n\n";
    
    std::cout << "Schwarzschild Radius: " << SCHWARZSCHILD_RADIUS << " units\n";
    std::cout << "Image Resolution: " << IMAGE_WIDTH << "x" << IMAGE_HEIGHT << "\n";
    std::cout << "Camera Distance: " << CAMERA_DISTANCE << " units\n\n";
    
    BlackHoleVisualizer visualizer;
    
    // Render the black hole
    visualizer.render();
    
    // Save the main visualization
    visualizer.saveToPPM("black_hole_visualization.ppm");
    
    // Generate spacetime curvature data
    visualizer.renderSpacetimeCurvature("spacetime_data");
    
    std::cout << "\nVisualization complete!\n";
    std::cout << "Files generated:\n";
    std::cout << "- black_hole_visualization.ppm (main image)\n";
    std::cout << "- spacetime_data_curvature.txt (curvature data)\n\n";
    
    std::cout << "To view the PPM image, use an image viewer that supports PPM format,\n";
    std::cout << "or convert it to PNG/JPEG using ImageMagick:\n";
    std::cout << "convert black_hole_visualization.ppm black_hole.png\n\n";
    
    std::cout << "The curvature data can be visualized with Python/matplotlib or gnuplot.\n";
    
    return 0;
}