# Tiny Learning-Based Model Predictive Control (Tiny-LB-MPC)

Tiny-LB-MPC is a novel framework for deploying learning-based model predictive control (MPC) on resource-constrained micro multirotor platforms. This repository contains the implementation used in the paper:

> **Tiny Learning-Based MPC for Multirotors: Solver-Aware Learning for Efficient Embedded Predictive Control**  
> *Babak Akbari, Justin Frank, Melissa Greeff*  
> [arXiv:2410.23634](https://arxiv.org/abs/2410.23634)

---

## Overview

Tiny aerial robots offer immense promise for applications like search and rescue, exploration, and monitoring. However, their computational constraints make it difficult to deploy advanced control algorithms.

Tiny-LB-MPC tackles this by:
- Leveraging the **differential flatness** property of multirotors to reduce the problem complexity.
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

### Clone the Repository
```bash
git clone https://github.com/BabakAkbari/Tiny-LB-MPC.git
```

### Prerequisites
- Arduino IDE or PlatformIO
- Teensyduino add-on
- CMake toolchain (optional)

### Build and Flash
- Open the firmware folder in the Arduino IDE.
- Select **Teensy 4.0** as the board.
- Compile and upload the code to your Teensy.

Refer to the `/firmware` and `/models` folders for code organization.

---

## Usage

- Load and initialize the MPC module.
- Send trajectory waypoints from the base station (or use predefined ones).
- Observe the performance through Crazyflie's logging and visualizations.

---

## Results

Tiny-LB-MPC achieves:
- **23% improvement in RMS tracking error** over a baseline MPC.
- Reliable **real-time operation at 100 Hz**.
- Full onboard computation without offloading.

---

## Citation

If you use Tiny-LB-MPC in your work, please cite:
```bibtex
@article{akbari2024tiny,
  title={Tiny Learning-Based MPC for Multirotors: Solver-Aware Learning for Efficient Embedded Predictive Control},
  author={Akbari, Babak and Frank, Justin and Greeff, Melissa},
  journal={arXiv preprint arXiv:2410.23634},
  year={2024}
}
```

---

## License

This repository is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

We thank the members of the Agile Robotics and Perception Lab for their support and feedback.
