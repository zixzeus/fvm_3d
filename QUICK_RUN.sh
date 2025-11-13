#!/bin/bash
# Quick script to build and run 3D Harris Sheet Magnetic Reconnection

set -e  # Exit on error

echo "╔════════════════════════════════════════════════════════════╗"
echo "║  FVM3D 3D Harris Sheet Magnetic Reconnection Simulator     ║"
echo "║  Quick Build & Run Script                                  ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# Configuration
BUILD_DIR="build"
EXECUTABLE="harris_sheet_3d"
NX=${1:-64}
NY=${2:-32}
NZ=${3:-32}

echo "═══════════════════════════════════════════════════════════════"
echo "Configuration:"
echo "  Build directory: ${BUILD_DIR}"
echo "  Resolution:     ${NX} × ${NY} × ${NZ}"
echo "═══════════════════════════════════════════════════════════════"
echo ""

# Step 1: Check if build directory exists
if [ ! -d "${BUILD_DIR}" ]; then
    echo "Creating build directory..."
    mkdir -p ${BUILD_DIR}
fi

# Step 2: Configure with CMake
echo "Configuring CMake..."
cd ${BUILD_DIR}
cmake -DCMAKE_BUILD_TYPE=Release ..
cd ..

# Step 3: Build
echo ""
echo "Building (this may take 1-2 minutes)..."
cd ${BUILD_DIR}
make -j4 ${EXECUTABLE}
cd ..

# Step 4: Run
echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║  Running 3D Magnetic Reconnection Simulation               ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

./${BUILD_DIR}/examples/${EXECUTABLE} ${NX} ${NY} ${NZ}

echo ""
echo "✓ Simulation complete!"
echo ""
echo "Results interpretation:"
echo "  max|B_y| - Should increase from 0.03 (reconnection indicator)"
echo "  div(B)   - Should stay < 1e-9 (GLM constraint maintained)"
echo "  KE, BE   - Energy redistribution during reconnection"
echo ""
echo "For more information, see RUN_3D_RECONNECTION.md"
