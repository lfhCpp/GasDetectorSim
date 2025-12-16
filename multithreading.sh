#!/bin/bash

# Check if parameters are correct
if [ $# -ne 1 ]; then
    echo "Usage: $0 <maximum build number N>"
    exit 1
fi

max_build=$1

# Clean up and create temporary directory
rm -rf tem
mkdir tem || exit 1
cd tem || exit 1

# Loop to create build directories
for i in $(seq 0 $max_build); do
    build_dir="build${i}"
    
    # Create build directory
    if ! mkdir "$build_dir"; then
        echo "Unable to create directory: $build_dir"
        continue
    fi

    # Execute build commands in a subshell
    (
        echo "Processing $build_dir..."
        cd "$build_dir" || exit 1
        
        # Execute cmake
        if ! cmake ../..; then
            echo "cmake failed: $build_dir"
            exit 1
        fi
        
        # Execute make
        if ! make -j8; then
            echo "make failed: $build_dir"
            exit 1
        fi
        
        echo "Starting MUSIC in background: $build_dir"
        ./MUSIC -m run.mac > "${build_dir}_music.log" 2>&1 &
        
        exit 0
    )
    
    # Check subshell execution status
    if [ $? -eq 0 ]; then
        echo "$build_dir processed successfully"
    else
        echo "$build_dir processing failed"
    fi
done

cd ..  # Return to initial directory

echo "All builds processed"