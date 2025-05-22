# Simulation of "A Tb/s Indoor MIMO Optical Wireless Backhaul System Using VCSEL Arrays"

This repository contains my personal simulation and implementation of the key ideas and system model described in the following research paper:

> *"A Tb/s indoor MIMO optical wireless backhaul system using VCSEL arrays."*  
> IEEE Transactions on Communications, vol. 70, no. 6, pp. 3995–4012, 2022.  
> [DOI: 10.1109/TCOMM.2022.3168035](https://doi.org/10.1109/TCOMM.2022.3168035)

---

## Purpose

This project aims to replicate and better understand the system performance and modeling techniques used in the above paper, particularly:

- High-speed **MIMO optical wireless communication**
- Use of **VCSEL arrays** for indoor backhaul
- Channel modeling, SINR calculations, and achievable data rate evaluation
- Spatial multiplexing and alignment strategies

> **Note:** This project is for **educational and research purposes only.** I am not affiliated with the original authors or their institutions.

---

## Reference

Kazemi, H., A. Khalighi, P. Leon, and S. Bourennane.  
*"A Tb/s indoor MIMO optical wireless backhaul system using VCSEL arrays."*  
IEEE Transactions on Communications 70.6 (2022): 3995–4012.

---

## Disclaimer

This is an independent, non-commercial simulation. The original paper and all credit belong to the authors. If you are interested in the original work, please refer to their publication on [IEEE Xplore](https://ieeexplore.ieee.org/document/9760970).

---

## Structure

```bash
/VCSEL_MIMO_Backhaul_Simulation
├── MIMO_backhaul.m
├── MIMO_displacement.m
├── PD_array.m
├── MIMO_FF.m
├── TEST.m
└── README.md
