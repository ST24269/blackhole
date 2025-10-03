#!/bin/bash

# Black Hole Visualization Build and Run Script
# This script compiles and runs the black hole visualization program

echo "========================================="
echo "Black Hole Spacetime Visualization"
echo "========================================="
echo

# Check if g++ is available
if ! command -v g++ &> /dev/null; then
    echo "Error: g++ compiler not found!"
    echo "Please install a C++ compiler:"
    echo "  Ubuntu/Debian: sudo apt install build-essential"
    echo "  CentOS/RHEL:   sudo yum install gcc-c++"
    echo "  macOS:         xcode-select --install"
    echo "  Windows:       Install MinGW or Visual Studio"
    exit 1
fi

# Check if source file exists
if [ ! -f "black_hole_viz.cpp" ]; then
    echo "Error: black_hole_viz.cpp not found!"
    echo "Please save the C++ source code as 'black_hole_viz.cpp'"
    exit 1
fi

echo "Step 1: Compiling the program..."
echo "Command: g++ -std=c++11 -O3 -march=native -funroll-loops -o black_hole_viz black_hole_viz.cpp -lm"
echo

# Compile with optimizations
if g++ -std=c++11 -O3 -march=native -funroll-loops -o black_hole_viz black_hole_viz.cpp -lm; then
    echo "✅ Compilation successful!"
    echo
else
    echo "❌ Compilation failed. Trying with basic flags..."
    if g++ -std=c++11 -O2 -o black_hole_viz black_hole_viz.cpp -lm; then
        echo "✅ Compilation successful with basic optimization!"
        echo
    else
        echo "❌ Compilation failed completely. Please check your compiler setup."
        exit 1
    fi
fi

echo "Step 2: Running the visualization..."
echo "This may take several minutes depending on your CPU..."
echo

# Create output directory
mkdir -p output
cd output

# Run the program
if ../black_hole_viz; then
    echo
    echo "✅ Visualization complete!"
    echo
    
    # Check output files
    if [ -f "black_hole_visualization.ppm" ]; then
        echo "Generated files:"
        echo "  📁 $(pwd)/black_hole_visualization.ppm"
        ls -lh black_hole_visualization.ppm
        
        # Try to convert to PNG if ImageMagick is available
        if command -v convert &> /dev/null; then
            echo
            echo "Step 3: Converting PPM to PNG..."
            if convert black_hole_visualization.ppm black_hole_visualization.png; then
                echo "✅ PNG conversion successful!"
                echo "  📁 $(pwd)/black_hole_visualization.png"
                ls -lh black_hole_visualization.png
            else
                echo "⚠️  PNG conversion failed, but PPM file is available"
            fi
        else
            echo
            echo "💡 Install ImageMagick to convert PPM to PNG:"
            echo "   Ubuntu/Debian: sudo apt install imagemagick"
            echo "   macOS:         brew install imagemagick"
            echo "   Or use any image viewer that supports PPM format"
        fi
    fi
    
    if [ -f "spacetime_data_curvature.txt" ]; then
        echo "  📁 $(pwd)/spacetime_data_curvature.txt"
        echo "     ($(wc -l < spacetime_data_curvature.txt) data points)"
    fi
    
    echo
    echo "Step 4: Python visualization (optional)..."
    
    # Check if Python and required packages are available
    if command -v python3 &> /dev/null; then
        echo "Checking Python packages..."
        
        # Test numpy and matplotlib
        if python3 -c "import numpy, matplotlib; print('✅ NumPy and Matplotlib available')" 2>/dev/null; then
            echo
            echo "You can now run the Python visualization script:"
            echo "  python3 ../curvature_visualization.py"
            echo
            echo "Or run it automatically? (y/n)"
            read -p "Run Python visualization now? " -n 1 -r
            echo
            
            if [[ $REPLY =~ ^[Yy]$ ]]; then
                if [ -f "../curvature_visualization.py" ]; then
                    echo "Running Python visualization..."
                    python3 ../curvature_visualization.py
                else
                    echo "❌ curvature_visualization.py not found in parent directory"
                fi
            fi
        else
            echo "⚠️  Python packages missing. Install with:"
            echo "   pip3 install numpy matplotlib scipy"
        fi
    else
        echo "⚠️  Python3 not found. Install Python to run additional visualizations."
    fi
    
else
    echo "❌ Program execution failed!"
    exit 1
fi

echo
echo "========================================="
echo "Summary"
echo "========================================="
echo "🎯 The black hole visualization demonstrates:"
echo "   • Spacetime curvature around a black hole"
echo "   • Gravitational lensing of background stars"
echo "   • Event horizon and photon sphere effects"
echo "   • Accretion disk with temperature gradients"
echo
echo "📁 Output files are in: $(pwd)"
echo
echo "🔬 Physics implemented:"
echo "   • Schwarzschild metric (general relativity)"
echo "   • Relativistic ray tracing"
echo "   • Gravitational light bending"
echo "   • Simplified Doppler effects"
echo
echo "📚 Educational value:"
echo "   • Visualizes Einstein's general relativity"
echo "   • Shows how massive objects warp spacetime"
echo "   • Demonstrates why black holes appear black"
echo "   • Illustrates gravitational lensing effects"
echo
echo "🎨 To view your results:"
echo "   • Open PPM files with GIMP, ImageMagick, or online viewers"
echo "   • Use the Python script for interactive 3D visualizations"
echo "   • Experiment with different parameters in the source code"
echo
echo "✨ Visualization complete! Enjoy exploring spacetime! ✨"