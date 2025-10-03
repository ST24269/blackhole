#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <SFML/Graphics.hpp>

// Define PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Vec3 {
    double x, y, z;
    
    Vec3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
    
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator*(double t) const { return Vec3(x * t, y * t, z * t); }
    Vec3 operator/(double t) const { return Vec3(x / t, y / t, z / t); }
    
    double dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    double length() const { return sqrt(x * x + y * y + z * z); }
    Vec3 normalize() const { double l = length(); return l > 0 ? *this / l : Vec3(); }
};

class BlackHole {
public:
    Vec3 position;
    double mass;
    double schwarzschild_radius;
    
    BlackHole(Vec3 pos, double m) : position(pos), mass(m) {
        schwarzschild_radius = 2.0 * mass;
    }
    
    Vec3 gravitational_field(const Vec3& point) const {
        Vec3 r = point - position;
        double dist = r.length();
        
        if (dist < schwarzschild_radius * 1.1) {
            return Vec3();
        }
        
        double force_magnitude = mass / (dist * dist * dist);
        double relativistic_correction = 1.0 + 3.0 * schwarzschild_radius / (2.0 * dist);
        force_magnitude *= relativistic_correction;
        
        return r * (-force_magnitude);
    }
    
    bool hits_event_horizon(const Vec3& point) const {
        return (point - position).length() <= schwarzschild_radius;
    }
};

class LightRay {
public:
    Vec3 position;
    Vec3 direction;
    std::vector<Vec3> path;
    bool escaped;
    bool absorbed;
    double energy;
    sf::Color color;
    
    LightRay(Vec3 pos, Vec3 dir, sf::Color c = sf::Color::Yellow) 
        : position(pos), direction(dir.normalize()), 
          escaped(false), absorbed(false), energy(1.0), color(c) {
        path.push_back(position);
    }
    
    void integrate_step(const BlackHole& bh, double dt) {
        if (absorbed || escaped) return;
        
        if (bh.hits_event_horizon(position)) {
            absorbed = true;
            color = sf::Color::Red;
            return;
        }
        
        if (position.length() > 50.0) {
            escaped = true;
            color = sf::Color::Green;
            return;
        }
        
        // Runge-Kutta 4th order integration
        Vec3 k1_pos = direction;
        Vec3 k1_dir = bh.gravitational_field(position);
        
        Vec3 k2_pos = direction + k1_dir * (dt / 2.0);
        Vec3 k2_dir = bh.gravitational_field(position + k1_pos * (dt / 2.0));
        
        Vec3 k3_pos = direction + k2_dir * (dt / 2.0);
        Vec3 k3_dir = bh.gravitational_field(position + k2_pos * (dt / 2.0));
        
        Vec3 k4_pos = direction + k3_dir * dt;
        Vec3 k4_dir = bh.gravitational_field(position + k3_pos * dt);
        
        position = position + (k1_pos + k2_pos * 2.0 + k3_pos * 2.0 + k4_pos) * (dt / 6.0);
        direction = direction + (k1_dir + k2_dir * 2.0 + k3_dir * 2.0 + k4_dir) * (dt / 6.0);
        direction = direction.normalize();
        
        double dist = (position - bh.position).length();
        if (dist > bh.schwarzschild_radius) {
            energy *= (1.0 - bh.schwarzschild_radius / (2.0 * dist));
        }
        
        // Update color based on energy (redshift effect)
        int red = std::min(255, (int)(255 * energy));
        int green = std::min(255, (int)(255 * energy * 0.8));
        int blue = std::min(255, (int)(255 * energy * 0.6));
        color = sf::Color(red, green, blue);
        
        path.push_back(position);
    }
};

class BlackHoleVisualization {
private:
    sf::RenderWindow window;
    BlackHole blackhole;
    std::vector<LightRay> light_rays;
    sf::Font font;
    
    // Camera and view parameters
    double zoom = 20.0;
    Vec3 camera_offset;
    bool simulation_running = false;
    bool simulation_paused = false;
    
    // UI elements
    sf::Text info_text;
    sf::Text controls_text;
    sf::CircleShape blackhole_shape;
    sf::CircleShape event_horizon_shape;
    
    // Simulation parameters
    double time_step = 0.01;
    int simulation_speed = 1;
    
public:
    BlackHoleVisualization() : window(sf::VideoMode(1200, 800), "Black Hole Light Ray Simulation"),
                               blackhole(Vec3(0, 0, 0), 1.0) {
        window.setFramerateLimit(60);
        
        // Try to load font (SFML includes default font)
        if (!font.loadFromFile("arial.ttf")) {
            // Use default font or handle error
            std::cout << "Warning: Could not load arial.ttf, using default font" << std::endl;
        }
        
        setup_ui();
        setup_blackhole_graphics();
    }
    
    void setup_ui() {
        info_text.setFont(font);
        info_text.setCharacterSize(14);
        info_text.setFillColor(sf::Color::White);
        info_text.setPosition(10, 10);
        
        controls_text.setFont(font);
        controls_text.setCharacterSize(12);
        controls_text.setFillColor(sf::Color::Cyan);
        controls_text.setPosition(10, window.getSize().y - 150);
        controls_text.setString(
            "Controls:\n"
            "1-5: Preset configurations\n"
            "SPACE: Start/Pause simulation\n"
            "R: Reset simulation\n"
            "C: Clear all rays\n"
            "Mouse: Click to add ray\n"
            "+/-: Adjust black hole mass\n"
            "WASD: Move camera\n"
            "Scroll: Zoom in/out"
        );
    }
    
    void setup_blackhole_graphics() {
        // Black hole (invisible - represents singularity)
        blackhole_shape.setRadius(2);
        blackhole_shape.setFillColor(sf::Color::Black);
        blackhole_shape.setOrigin(2, 2);
        
        // Event horizon
        double horizon_radius = blackhole.schwarzschild_radius * zoom;
        event_horizon_shape.setRadius(horizon_radius);
        event_horizon_shape.setFillColor(sf::Color(50, 0, 0, 100));
        event_horizon_shape.setOutlineColor(sf::Color::Red);
        event_horizon_shape.setOutlineThickness(2);
        event_horizon_shape.setOrigin(horizon_radius, horizon_radius);
    }
    
    Vec3 screen_to_world(sf::Vector2i screen_pos) {
        sf::Vector2f world_pos = window.mapPixelToCoords(screen_pos);
        return Vec3((world_pos.x - window.getSize().x/2) / zoom + camera_offset.x,
                   (world_pos.y - window.getSize().y/2) / zoom + camera_offset.y, 0);
    }
    
    sf::Vector2f world_to_screen(const Vec3& world_pos) {
        double screen_x = (world_pos.x - camera_offset.x) * zoom + window.getSize().x/2;
        double screen_y = (world_pos.y - camera_offset.y) * zoom + window.getSize().y/2;
        return sf::Vector2f(screen_x, screen_y);
    }
    
    void handle_input() {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
            
            if (event.type == sf::Event::KeyPressed) {
                switch (event.key.code) {
                    case sf::Keyboard::Space:
                        if (simulation_running) {
                            simulation_paused = !simulation_paused;
                        } else {
                            simulation_running = true;
                            simulation_paused = false;
                        }
                        break;
                        
                    case sf::Keyboard::R:
                        reset_simulation();
                        break;
                        
                    case sf::Keyboard::C:
                        light_rays.clear();
                        break;
                        
                    case sf::Keyboard::Num1:
                        create_light_ring(8.0, 16);
                        break;
                        
                    case sf::Keyboard::Num2:
                        create_parallel_beam(15.0, 10, 5.0);
                        break;
                        
                    case sf::Keyboard::Num3:
                        create_spiral_pattern();
                        break;
                        
                    case sf::Keyboard::Num4:
                        create_cross_pattern();
                        break;
                        
                    case sf::Keyboard::Num5:
                        create_random_rays();
                        break;
                        
                    case sf::Keyboard::Equal:
                        adjust_black_hole_mass(1.2);
                        break;
                        
                    case sf::Keyboard::Hyphen:
                        adjust_black_hole_mass(0.8);
                        break;
                        
                    case sf::Keyboard::W:
                        camera_offset.y -= 1.0 / zoom;
                        break;
                        
                    case sf::Keyboard::S:
                        camera_offset.y += 1.0 / zoom;
                        break;
                        
                    case sf::Keyboard::A:
                        camera_offset.x -= 1.0 / zoom;
                        break;
                        
                    case sf::Keyboard::D:
                        camera_offset.x += 1.0 / zoom;
                        break;
                }
            }
            
            if (event.type == sf::Event::MouseButtonPressed) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    add_ray_from_mouse(sf::Vector2i(event.mouseButton.x, event.mouseButton.y));
                }
            }
            
            if (event.type == sf::Event::MouseWheelScrolled) {
                if (event.mouseWheelScroll.delta > 0) {
                    zoom *= 1.1;
                } else {
                    zoom /= 1.1;
                }
                zoom = std::max(1.0, std::min(100.0, zoom));
                setup_blackhole_graphics();
            }
        }
    }
    
    void add_ray_from_mouse(sf::Vector2i mouse_pos) {
        Vec3 start_pos = screen_to_world(mouse_pos);
        Vec3 direction = (Vec3(0, 0, 0) - start_pos).normalize();
        
        // Add some randomness to direction
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> noise(0.0, 0.1);
        direction = direction + Vec3(noise(gen), noise(gen), 0);
        direction = direction.normalize();
        
        light_rays.emplace_back(start_pos, direction, sf::Color::Cyan);
    }
    
    void create_light_ring(double radius, int num_rays) {
        light_rays.clear();
        for (int i = 0; i < num_rays; ++i) {
            double theta = 2.0 * M_PI * i / num_rays;
            Vec3 start_pos(radius * cos(theta), radius * sin(theta), 0);
            Vec3 direction = (Vec3(0, 0, 0) - start_pos).normalize();
            
            sf::Color ray_color = sf::Color(
                128 + 127 * cos(theta),
                128 + 127 * sin(theta),
                200
            );
            
            light_rays.emplace_back(start_pos, direction, ray_color);
        }
    }
    
    void create_parallel_beam(double distance, int num_rays, double beam_width) {
        light_rays.clear();
        for (int i = 0; i < num_rays; ++i) {
            double y = (i - num_rays/2.0) * beam_width / num_rays;
            Vec3 start_pos(-distance, y, 0);
            Vec3 direction(1, 0, 0);
            light_rays.emplace_back(start_pos, direction, sf::Color::Yellow);
        }
    }
    
    void create_spiral_pattern() {
        light_rays.clear();
        int num_rays = 20;
        for (int i = 0; i < num_rays; ++i) {
            double t = i * 0.5;
            double radius = 5 + t;
            double theta = t * 2;
            
            Vec3 start_pos(radius * cos(theta), radius * sin(theta), 0);
            Vec3 direction = Vec3(-sin(theta), cos(theta), 0);
            light_rays.emplace_back(start_pos, direction, sf::Color::Magenta);
        }
    }
    
    void create_cross_pattern() {
        light_rays.clear();
        int rays_per_arm = 8;
        double distance = 10;
        
        // Horizontal rays
        for (int i = 0; i < rays_per_arm; ++i) {
            double y = (i - rays_per_arm/2.0) * 0.5;
            light_rays.emplace_back(Vec3(-distance, y, 0), Vec3(1, 0, 0), sf::Color::Red);
            light_rays.emplace_back(Vec3(distance, y, 0), Vec3(-1, 0, 0), sf::Color::Red);
        }
        
        // Vertical rays
        for (int i = 0; i < rays_per_arm; ++i) {
            double x = (i - rays_per_arm/2.0) * 0.5;
            light_rays.emplace_back(Vec3(x, -distance, 0), Vec3(0, 1, 0), sf::Color::Blue);
            light_rays.emplace_back(Vec3(x, distance, 0), Vec3(0, -1, 0), sf::Color::Blue);
        }
    }
    
    void create_random_rays() {
        light_rays.clear();
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> angle_dist(0, 2 * M_PI);
        std::uniform_real_distribution<double> radius_dist(8, 15);
        
        for (int i = 0; i < 25; ++i) {
            double theta = angle_dist(gen);
            double radius = radius_dist(gen);
            
            Vec3 start_pos(radius * cos(theta), radius * sin(theta), 0);
            Vec3 direction = (Vec3(0, 0, 0) - start_pos).normalize();
            
            sf::Color color(rand() % 256, rand() % 256, rand() % 256);
            light_rays.emplace_back(start_pos, direction, color);
        }
    }
    
    void adjust_black_hole_mass(double factor) {
        double new_mass = blackhole.mass * factor;
        new_mass = std::max(0.1, std::min(10.0, new_mass));
        blackhole = BlackHole(Vec3(0, 0, 0), new_mass);
        setup_blackhole_graphics();
    }
    
    void reset_simulation() {
        simulation_running = false;
        simulation_paused = false;
        for (auto& ray : light_rays) {
            ray.absorbed = false;
            ray.escaped = false;
            ray.energy = 1.0;
            ray.position = ray.path.empty() ? Vec3() : ray.path[0];
            ray.path.clear();
            ray.path.push_back(ray.position);
            ray.color = sf::Color::Yellow;
        }
    }
    
    void update_simulation() {
        if (!simulation_running || simulation_paused) return;
        
        for (int step = 0; step < simulation_speed; ++step) {
            bool any_active = false;
            for (auto& ray : light_rays) {
                if (!ray.absorbed && !ray.escaped) {
                    ray.integrate_step(blackhole, time_step);
                    any_active = true;
                }
            }
            if (!any_active) {
                simulation_paused = true;
                break;
            }
        }
    }
    
    void render() {
        window.clear(sf::Color::Black);
        
        // Draw grid
        draw_grid();
        
        // Draw event horizon
        sf::Vector2f bh_screen = world_to_screen(blackhole.position);
        event_horizon_shape.setPosition(bh_screen);
        window.draw(event_horizon_shape);
        
        // Draw black hole center
        blackhole_shape.setPosition(bh_screen);
        window.draw(blackhole_shape);
        
        // Draw light rays
        for (const auto& ray : light_rays) {
            if (ray.path.size() > 1) {
                for (size_t i = 1; i < ray.path.size(); ++i) {
                    sf::Vector2f start = world_to_screen(ray.path[i-1]);
                    sf::Vector2f end = world_to_screen(ray.path[i]);
                    
                    sf::Vertex line[] = {
                        sf::Vertex(start, ray.color),
                        sf::Vertex(end, ray.color)
                    };
                    window.draw(line, 2, sf::Lines);
                }
            }
        }
        
        // Update and draw UI
        update_info_text();
        window.draw(info_text);
        window.draw(controls_text);
        
        window.display();
    }
    
    void draw_grid() {
        sf::Color grid_color(50, 50, 50, 128);
        
        // Draw grid lines
        for (double x = -50; x <= 50; x += 5) {
            sf::Vector2f start = world_to_screen(Vec3(x, -50, 0));
            sf::Vector2f end = world_to_screen(Vec3(x, 50, 0));
            
            if (start.x >= 0 && start.x <= window.getSize().x) {
                sf::Vertex line[] = {
                    sf::Vertex(start, grid_color),
                    sf::Vertex(end, grid_color)
                };
                window.draw(line, 2, sf::Lines);
            }
        }
        
        for (double y = -50; y <= 50; y += 5) {
            sf::Vector2f start = world_to_screen(Vec3(-50, y, 0));
            sf::Vector2f end = world_to_screen(Vec3(50, y, 0));
            
            if (start.y >= 0 && start.y <= window.getSize().y) {
                sf::Vertex line[] = {
                    sf::Vertex(start, grid_color),
                    sf::Vertex(end, grid_color)
                };
                window.draw(line, 2, sf::Lines);
            }
        }
    }
    
    void update_info_text() {
        int active = 0, absorbed = 0, escaped = 0;
        for (const auto& ray : light_rays) {
            if (ray.absorbed) absorbed++;
            else if (ray.escaped) escaped++;
            else active++;
        }
        
        std::stringstream ss;
        ss << std::fixed << std::setprecision(2);
        ss << "Black Hole Mass: " << blackhole.mass << "\n";
        ss << "Schwarzschild Radius: " << blackhole.schwarzschild_radius << "\n";
        ss << "Zoom: " << zoom << "x\n";
        ss << "Light Rays: " << light_rays.size() << "\n";
        ss << "Active: " << active << " | Absorbed: " << absorbed << " | Escaped: " << escaped << "\n";
        ss << "Status: " << (simulation_running ? (simulation_paused ? "PAUSED" : "RUNNING") : "STOPPED");
        
        info_text.setString(ss.str());
    }
    
    void run() {
        std::cout << "Black Hole Simulation Started!" << std::endl;
        std::cout << "Use number keys 1-5 for preset configurations" << std::endl;
        
        while (window.isOpen()) {
            handle_input();
            update_simulation();
            render();
        }
    }
};

int main() {
    std::cout << "=== BLACK HOLE LIGHT RAY SIMULATION (GUI) ===" << std::endl;
    std::cout << "Starting graphical simulation..." << std::endl;
    
    try {
        BlackHoleVisualization sim;
        sim.run();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}