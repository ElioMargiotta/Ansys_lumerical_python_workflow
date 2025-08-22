import subprocess
import time

def run_script(path):
    start_time = time.perf_counter()
    print(f"\n🟩 Running: {path}")
    result = subprocess.run(["python", path], capture_output=True, text=True)
    end_time = time.perf_counter()

    elapsed = end_time - start_time
    if result.returncode == 0:
        print(f"✅ Done: {path} in {elapsed:.2f} seconds")
        print(result.stdout)
    else:
        print(f"❌ Error in {path}:")
        print(result.stderr)
        exit(1)

if __name__ == "__main__":
    print("🧪 Starting Optical Simulation Pipeline")
    total_start = time.perf_counter()

    # Step 1: Build structure and run simulation
    run_script("scripts/build_structure.py")

    # Step 2: Analyze results
    run_script("scripts/analyse_results.py")

    total_end = time.perf_counter()
    total_elapsed = total_end - total_start

    print(f"\n🎉 All steps completed successfully in {total_elapsed:.2f} seconds.")
