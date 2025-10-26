# NiTi Nanoparticle SME/SE Simulation Workspace

[](https://github.com/PrudhviNallagatla/SME_SE)
[](https://www.lammps.org)
[](https://www.python.org/downloads/)

This repository provides a complete, scalable workflow for generating and simulating **Nickel-Titanium (NiTi)** nanoparticles to study the **Shape Memory Effect (SME)** and **Superelasticity (SE)** using LAMMPS.

The primary focus is an optimized workflow for parallel simulations on **Azure Virtual Machines**.

## Features

  * **Nanoparticle Generation**: Create realistic NiTi nanoparticles with various shapes (sphere, blob, ellipsoid, faceted, rough), structures (single crystal, polycrystalline, amorphous), and defects (vacancies, antisites, surface vacancies).
  * **Simulation Pipelines**: Automated scripts for energy minimization, thermal cycling, SME deformation, and SE mechanical loading.
  * **Analysis Tools**: Jupyter notebooks and Python scripts for analyzing transformation temperatures ($M_s, M_f, A_s, A_f$) and mechanical properties from simulation data.
  * **EDM Synthesis Focus**: Optimized for simulating Electrical Discharge Machining (EDM) conditions with rapid quench and amorphous structures.

-----

## Navigating This Repository

This repository is best navigated using a **tree view** (e.g., in VS Code's file explorer) to understand the directory structure.

All primary simulation workflows are controlled by shell scripts located in the [`/scripts`](./scripts/) directory. The core LAMMPS input files are in [`/src/lmps`](./src/lmps/), and the Python-based nanoparticle generator is [`/src/niti_np_gen.py`](./src/niti_np_gen.py).

## File Structure

```
.
├── all_dependencies.sh           # Installation script for LAMMPS and Python deps
├── README.md                     # This file
├── data/
│   ├── analysis_output/          # Analysis results (plots, summaries)
│   ├── raw_output/               # Simulation outputs (minimize, thermal_cycle, sme, se)
│   └── structures/               # Generated nanoparticle data files
├── logs/                         # Simulation log files
├── potentials/                   # MEAM potential files (NiTi.meam, library.meam)
├── scripts/                      # *** ALL WORKFLOW SCRIPTS ARE HERE ***
│   ├── run_all.sh                # Full pipeline script (SME + SE)
│   ├── structure_gen.sh          # Structure generation wrapper
│   ├── se_pipeline.sh            # SE simulation pipeline
│   └── sme_scripts/
│       ├── run_01_characterize.sh  # SME (minimize + thermal cycle)
│       ├── run_02_deform.sh        # SME (deformation)
│       └── analyze_thermal_cycle.ipynb  # Analysis notebook
├── src/
│   ├── niti_np_gen.py            # Nanoparticle generator (Python)
│   └── lmps/                     # LAMMPS input scripts
└── .devcontainer/                # Docker config for LOCAL testing
```

-----

## Recommended Azure Workflow Overview

This is our primary, scalable workflow. It's designed to run many simulations in parallel cheaply using Azure Spot Instances.

```
 [ 1. Build "Golden Image" (One-Time) ]
 (Ubuntu + LAMMPS + Python + Ovito)
                |
                v
 [ 2. Set Up Azure File Share (One-Time) ]
 (Upload /src, /potentials, /scripts, /data)
                |
                v
 [ 3. Launch N VMs from "Golden Image" ]
  /                  |                \
 /                   |                 \
v                    v                  v
[ VM 1 ]            [ VM 2 ]            [ VM 3 ]
(Mount Share)       (Mount Share)       (Mount Share)
(Run Jobs 1-4)      (Run Jobs 5-8)      (Run Jobs 9-10)
  |                    |                  |
  \                    v                  /
   \          (Write Results)            /
    \                 |                 /
     `------> [ /workspace (Shared) ] <----'
                (Azure File Share)
                |
                v
 [ 4. Analyze Results (on one VM) ]
 (Jupyter + VS Code)
```

-----

## Detailed Storage Setup 🗂️

Imagine on your **local computer** (your laptop), your project is in a folder named `NiTi_Project`. When you open it, you see this structure:

```
NiTi_Project/
├── .git/
├── all_dependencies.sh
├── README.md
├── data/
│   ├── structures/
│   └── ...
├── logs/
├── potentials/
│   ├── NiTi.meam
│   └── library.meam
├── scripts/
│   ├── run_all.sh
│   └── ...
└── src/
    ├── niti_np_gen.py
    └── ...
```

Your goal is to make the *inside* of your Azure File Share look exactly like the *inside* of this `NiTi_Project` folder.

### Step 1: Go to Your File Share in Azure

1.  In the Azure Portal, find and click on your **Storage Account**.
2.  On the left-hand menu, click on **"File shares"**.
3.  Click on the name of the file share you created (e.g., `rimuru-workspace-storage`).

You are now looking at the **root** of your file share. It is currently empty. You will see an **"Upload"** button at the top.

### Step 2: Upload Your Project Folders

1.  Click the **"Upload"** button. A new panel will open.
2.  On your **local computer**, open your `NiTi_Project` folder.
3.  **Select** the following five folders (and *only* these five folders):
      * `data`
      * `logs`
      * `potentials`
      * `scripts`
      * `src`
4.  **Drag and drop** these five folders from your computer into the Azure "Upload" panel.

**Note:** You do not need to upload `all_dependencies.sh`, `README.md`, or the `.git` folder. Your "Template Image" *already* has all the dependencies installed, and the simulation scripts don't need the other files.

### Step 3: Verify the Final Structure

After the upload finishes, the root directory of your Azure File Share should look like this:

```
<Your File Share Name> /
├── data/
├── logs/
├── potentials/
├── scripts/
└── src/
```

When you (for example) click on the `potentials` folder *in the Azure Portal*, you should see your `NiTi.meam` and `library.meam` files inside it.

### Why This Works

Now, when you do the **simulation step**, you will run this command on your new VMs:

```bash
# 2a. Create the empty "link" folder
mkdir /home/rimuru/workspace

# 2b. Mount your storage to that folder
sudo mount -t cifs //yourstorage.file.core.windows.net/... /home/rimuru/workspace -o <...options...>
```

You have now told the VM: "Whenever any program (like LAMMPS or Python) tries to access the `/home/rimuru/workspace` directory, you must *actually* go to the Azure File Share."

So, when your script (which you never had to edit!) runs a command:
`lmp -in /home/rimuru/workspace/src/lmps/sme_thermal_cycle.lmp ...`

The VM does the following:

1.  It goes to `/home/rimuru/workspace`...
2.  ...realizes this is just a link to your **Azure File Share**...
3.  ...opens the `src` folder on the file share...
4.  ...opens the `lmps` folder...
5.  ...and finds your `sme_thermal_cycle.lmp` file.

When LAMMPS writes an output file to `/home/rimuru/workspace/data/raw_output/`, it is writing *directly* onto your permanent Azure File Share, not the temporary VM.


* Which is faster, 8x2 or 4x4, is a classic performance-tuning question. There is no single "right" answer; it depends on your exact simulation:

    **8x2 (More MPI)**: Often better for very large systems where splitting the domain into more pieces is most important.

    **4x4 (More OMP)**: Often better for systems with complex calculations (like MEAM or REAXFF), where giving more thread-power to each MPI rank to calculate forces is more important.

-----

## Usage

### 1\. Local / Test Setup (Docker)

The included Docker configuration (`.devcontainer/`) is intended for **local testing and development** on low-end PCs or for use on HPCs where you cannot create custom VM images. It is *not* the recommended workflow for large-scale runs.

1.  **Clone the workspace** and open it in a dev container (e.g., in VS Code with the "Dev Containers" extension).
2.  **Install dependencies** (this will run automatically in the dev container):
    ```bash
    ./all_dependencies.sh
    ```
3.  **Verify installation**:
      * `lmp -help`
      * `python3 --version`

### 2\. Cloud Simulation Workflow (Azure)

This is the most powerful and cost-effective method for running your simulations.

#### Step 1: Build Your "Golden" VM Image (One-Time Setup)

This step saves you hours of setup time for every run. You do this **once**.

1.  **Launch the Base VM:**
      * **Image:** Ubuntu Server 22.04 LTS
      * **Size:** `Standard_F16s_v2` (16 cores, 32GB RAM)
      * **Pricing:** **Azure Spot Instance** (set eviction to "Deallocate")
      * **Username:** e.g., `rimuru`
      * **OS Disk:** **32 GB Premium SSD (P4)**
2.  **Connect & Install:**
      * Connect to the VM using **VS Code Remote-SSH** (`ssh rimuru@Your_VM_IP`).
      * Upload or paste your `all_dependencies.sh` script.
      * Run the installation:
        ```bash
        chmod +x all_dependencies.sh
        ./all_dependencies.sh
        ```
      * This will take **20-40 minutes** to compile LAMMPS and install all Python packages. ☕
3.  **Clean Up & Shut Down:**
      * Run `history -c` (optional).
      * Shut down the VM cleanly: `sudo shutdown now`
4.  **Capture the Image:**
      * In the Azure Portal, find your (now stopped) VM.
      * Click the **"Capture"** button.
      * Give the image a name (e.g., `rimuru-lammps-ubuntu-v1`). This is now your reusable "Golden Image."

#### Step 2: Set Up Shared Storage (One-Time Setup)

1.  **Create Storage:** In the Azure Portal, create a **256 GB Azure File Share**.
2.  **Upload Project:** Upload your *entire* project folder to this file share. It should now contain `/data`, `/potentials`, `/src`, and `/scripts`.

#### Step 3: Run Simulations in Parallel 🚀

### Before starting, very important: VM Management and Cost Optimization 💰

**🚀 Flexible VM Scheduling:** It doesn't matter whether you launch all VMs simultaneously or run them sequentially (one after another) - Azure charges per VM per hour, so the total cost remains the same. This gives you maximum flexibility: run all at once for fastest results, or stagger them to match your schedule or budget constraints.

**⚡ Golden Image Priority:** Always set up your Golden VM Image correctly **first** before starting any simulations. Once ready, you can "whip out" VMs instantly - launch, run jobs, and delete. Each simulation cycle is: VM up → mount storage → run → VM down. No setup time wasted!

**🧠 Don't Overpay for RAM: Core-First VM Selection:** **Cores are your top priority** for parallel LAMMPS simulations. For less demanding jobs, choose VMs with more cores but less RAM (e.g., 32 cores with 16GB RAM) over fewer cores with excessive RAM. Your Golden Image boots in minutes, so focus on computational power - permanent storage handles all data persistence.

**💸 Late-Night Optimization:** If running simulations late at night, remember: only storage costs matter. VM compute costs stop when you delete the VMs. Keep your Azure File Share for long-term data, but delete VMs immediately after jobs complete to minimize expenses.

**VM Selection Tip (very important):** When selecting VMs for simulations, prioritize configurations with higher core counts over excessive RAM, as permanent storage handles data persistence. A golden VM image ensures quick setup, so focus on cores for parallel processing—32GB RAM is typically sufficient, and even 16GB may suffice in some cases, but never compromise on core availability.
      * These VMs will boot in minutes, fully configured.
  
Now, whenever you need to run your study:

1.  **Launch VMs:** Launch **3** (or more) new VMs.
      * **Image:** Under "Image," select **"My Images"** and choose your `rimuru-lammps-ubuntu-v1` image.
      * **Settings:** Use the same `Standard_F16s_v2`, **Spot Instance**, and **32GB P4** disk settings.
2.  **Connect & Mount:**
      * Open a separate VS Code window (Remote-SSH) for each VM.
      * In *each* VM's terminal, mount your shared storage (Azure provides the copy-paste command for this):
        ```bash
        sudo mkdir /workspace
        sudo mount -t cifs //yourstorage.file.core.windows.net/... /workspace
        ```
3.  **Assign Jobs:**
      * **On VM 1:**
          * Copy the pipeline script locally: `cp /workspace/scripts/run_all.sh ~/run_vm1.sh`
          * Edit `~/run_vm1.sh` in VS Code to run *only* the 5nm and 10nm files.
          * Run it: `bash ~/run_vm1.sh`
      * **On VM 2:**
          * Copy the script: `cp /workspace/scripts/run_all.sh ~/run_vm2.sh`
          * Edit `~/run_vm2.sh` to run *only* the 15nm and 20nm files.
          * Run it: `bash ~/run_vm2.sh`
      * **On VM 3:**
          * Copy the script: `cp /workspace/scripts/run_all.sh ~/run_vm3.sh`
          * Edit `~/run_vm3.sh` to run *only* the 25nm files.
          * Run it: `bash ~/run_vm3.sh`
4.  **Wait & Shut Down:**
      * All VMs are now running in parallel, reading from and writing to your central `/workspace` storage.
      * When the jobs are finished, **Delete the VMs** to stop all costs. Your results are safe on the shared storage.

### 3\. Key Scripts and Usage

#### Generate Structures

Use [`src/niti_np_gen.py`](src/niti_np_gen.py) to generate your nanoparticles *before* running simulations.

```bash
# Generate a 10nm amorphous blob and its LAMMPS data file
python3 src/niti_np_gen.py --diameter 10 --shape blob --amorphous 0.3 --output niti_10nm_edm.xyz --lammps --seed 42

# Generate a 15nm single crystal sphere
python3 src/niti_np_gen.py --diameter 15 --shape sphere --output niti_15nm_crystal.xyz --lammps
```

#### Run Full Pipeline (SME + SE)

The [`scripts/run_all.sh`](scripts/run_all.sh) script is a master script that processes all `.data` files it finds in [`data/structures`](./data/structures/).

```bash
# This will find all structures and run the full SME and SE pipelines on them
./scripts/run_all.sh
```

*As shown in the Azure workflow, you should **copy and edit** this script on each VM to divide the workload.*

#### Run Individual Pipelines

  * **SME Pipeline:**
    1.  [`scripts/sme_scripts/run_01_characterize.sh`](scripts/sme_scripts/run_01_characterize.sh) (Minimizes + runs thermal cycle)
    2.  [`scripts/sme_scripts/run_02_deform.sh`](scripts/sme_scripts/run_02_deform.sh) (Deforms the particle)
  * **SE Pipeline:**
      * [`scripts/se_pipeline.sh`](scripts/se_pipeline.sh) (Minimizes + runs mechanical loading/unloading)

-----

### ⚠️ Important: Data Reduction for Storage

For this type of study, you **do not** need multi-gigabyte dump files. The primary results come from small `thermo_detailed.txt` files.

To save space (and stay under 256 GB), you **must** edit your LAMMPS scripts in [`src/lmps/`](./src/lmps/) (e.g., `sme_thermal_cycle.lmp`, `se_load.lmp`).

Find the `dump` commands and increase the output interval to a very large number:

  * **FROM:** `dump 1 all custom 2000 ${output_dir}/sme_*.dump ...`
  * **TO:** `dump 1 all custom 500000 ${output_dir}/sme_*.dump ...`

This saves only a few snapshots for figures, not thousands of frames for a video.

### 💡 Pro-Tip for Cost Savings while zipping

You **do not need an expensive 16-core VM** to run the compression script.

Since the task is limited by storage I/O, a powerful CPU will just sit idle. You can save money by launching a cheap 2-core or 4-core VM just for this task. It will likely take the exact same amount of time (30-45 minutes) and cost you pennies.
-----

## 4\. Data Analysis Workflow 📊

You can analyze all your data in the cloud using VS Code's remote Jupyter capabilities.

1.  **Keep One VM Running** (or launch a new one from your Golden Image).
2.  **Connect** via VS Code Remote-SSH and **mount** your `/workspace` shared drive.
3.  **Launch Jupyter:** In the VS Code terminal:
    ```bash
    jupyter notebook --no-browser --port=8888
    ```
4.  **Connect VS Code to Jupyter:**
      * Press `F1` and type `Jupyter: Specify Jupyter Server for Connections`.
      * Select "Existing" and paste the `http://localhost:8888/?token=...` URL from your terminal.
5.  **Analyze\!**
      * Open [`scripts/sme_scripts/analyze_thermal_cycle.ipynb`](./scripts/sme_scripts/analyze_thermal_cycle.ipynb).
      * Select your remote server as the kernel (top-right corner).
      * You can now run Python code *on the VM* with direct access to all your results in `/workspace/data/raw_output/`.

### Analysis Objectives

  * **SME Analysis:**
      * Plot **Box Length vs. Temperature** from the thermal cycling data to visualize the SME hysteresis loop and recovery.
      * Calculate **Shape Recovery % vs. Particle Size** to quantify the size-dependent effect.
  * **SE Analysis:**
      * Plot **Stress vs. Strain** from the mechanical loading data to observe the characteristic superelastic "flag" plateau.
  * **General:**
      * Plot **Potential Energy vs. Temperature** to identify transformation temperatures ($M_s, M_f, A_s, A_f$) using the derivative method.
      * Use **Ovito** (which is installed by `all_dependencies.sh`) on your few saved dump files to create images showing the Austenite (B2) and Martensite (B19') phases.

## Dependencies

  * **LAMMPS**: Compiled with MEAM, SNAP, and other packages (installed via `all_dependencies.sh`)
  * **Python 3**: With NumPy, SciPy, Matplotlib, Pandas, Jupyter, Ovito
  * **System**: Ubuntu 22.04, MPI (mpirun), Git

## Troubleshooting

  * **LAMMPS errors**: Check that the potential files exist in [`/potentials`](./potentials) and are correctly referenced.
  * **Analysis failures**: You may need to adjust `ANALYSIS_CONFIG` parameters in the analysis notebooks (e.g., increase `WINDOW_LENGTH` for noisy data).
  * **Azure Spot Eviction**: If a Spot VM is evicted (deallocated), your job will stop. Simply re-launch a new VM from your Golden Image, mount the storage, and re-run the script for the *missing* data. No results will be lost.