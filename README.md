# Tiny Learning-Based Model Predictive Control (Tiny-LB-MPC)

Tiny-LB-MPC is a novel framework for deploying learning-based model predictive control (MPC) on resource-constrained micro multirotor platforms.


## Overview

Tiny-LB-MPC tackles this by:
- Leveraging the **differential flatness** property of multirotors to reduce problem complexity.
- Introducing a **solver-aware learning** approach that makes the control problem tractable on embedded hardware.
- Achieving **real-time control at 100 Hz** onboard microcontrollers.

This is the **first implementation of learning-based MPC deployed onboard** a tiny (53 g) multirotor system.

## Hardware Setup

- **Crazyflie 2.1**
- **Teensy 4.0** (executes MPC)
- **Custom expansion board** (to interface Crazyflie and Teensy)
- **High-thrust motors and propellers**

## Getting Started

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
- **43% improvement in RMS tracking error** over baseline MPC.
Here are three videos comparing Tiny LB MPC with Tiny L MPC and Tiny FB MPC:

<p align="center">
  <img src="media/Tiny L MPC.gif" width="300" alt="Tiny L MPC">
  <img src="media/Tiny FB MPC.gif" width="300" alt="Tiny FB MPC">
  <img src="media/Tiny LB MPC.gif" width="300" alt="Tiny LB MPC">
</p>

### Citation

If you use Tiny-LB-MPC in your work, please cite:

TBD
