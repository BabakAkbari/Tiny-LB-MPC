# Tiny Learning-Based Model Predictive Control (Tiny-LB-MPC)

Tiny-LB-MPC is a novel framework for deploying learning-based model predictive control (MPC) on resource-constrained micro multirotor platforms.

---

## Overview

Tiny aerial robots offer immense promise for applications like search and rescue, exploration, and monitoring. However, their computational constraints make it difficult to deploy advanced control algorithms.

Tiny-LB-MPC tackles this by:
- Leveraging the **differential flatness** property of multirotors to reduce problem complexity.
- Introducing a **solver-aware learning** approach that makes the control problem tractable on embedded hardware.
- Achieving **real-time control at 100 Hz** onboard microcontrollers.

This is the **first implementation of learning-based MPC deployed onboard** a tiny (53 g) multirotor system.

---

## Features

- **Learning-Based MPC**: Augments model-based control with data-driven dynamics.
- **Solver-Aware Learning**: Optimizes for runtime efficiency during training.
- **Embedded Deployment**: Runs entirely on a Teensy 4.0 microcontroller.
- **Real-Time Capable**: Executes at 100 Hz using limited onboard resources.
- **Tracking Improvement**: Demonstrates ~23% improvement over baseline MPC methods.

---

## Hardware Setup

- **Crazyflie 2.1** (flight controller and sensors)
- **Teensy 4.0** (executes MPC and low-level control)
- **Custom expansion board** (to interface Crazyflie and Teensy)
- **High-thrust motors and propellers** (to support aggressive flight)

---

## Getting Started

TBD

### Prerequisites

- [Visual Studio Code](https://code.visualstudio.com/)
- [PlatformIO VS Code extension](https://platformio.org/install/ide?install=vscode)
- Teensy platform for PlatformIO

### Teensy Setup for PlatformIO

1. Install PlatformIO in VS Code.
2. Install the Teensy platform:
   ```bash
   pio platform install teensy
   ```

## Results

Tiny-LB-MPC achieves:
- **23% improvement in RMS tracking error** over baseline MPC.
- **Reliable real-time operation at 100 Hz**.
- Fully embedded deployment—no offboard computation required.

### GIFs

Here are three GIFs showcasing different versions of the Tiny-LB-MPC implementation:

<div>
  <img src="media/Tiny L MPC.gif" width="290" alt="Tiny L MPC" style="display: inline-block;">
  <img src="media/Tiny FB MPC.gif" width="290" alt="Tiny FB MPC" style="display: inline-block;">
  <img src="media/Tiny LB MPC.gif" width="290" alt="Tiny LB MPC" style="display: inline-block;">
</div>

## Citation

If you use Tiny-LB-MPC in your work, please cite:

TBD

---

## License

This repository is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
