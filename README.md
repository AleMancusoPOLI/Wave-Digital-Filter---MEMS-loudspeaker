# Wave Digital Filter (WDF) Modeling of a Piezoelectric MEMS Loudspeaker

## Assignment Overview 
In this assignment, we developed a Wave Digital Filter (**WDF**) model of a piezoelectric **MEMS** loudspeaker, based on its linear lumped-element model (**LEM**) representation.
Starting from the equivalent mechanical domain circuit, the corresponding WDF was derived according to the following schematic:

![Immagine WhatsApp 2025-06-09 ore 22 43 22_6c010247](https://github.com/user-attachments/assets/4cda7e78-a396-40de-8293-9a86d7285b91)

## Repository Contents
| File                   | Description                                           |
|-----------------------|-------------------------------------------------------|
| `mems_spk_WD.m`          | The main MATLAB script implementing the WDF model.    |
| `mems_spk_SSC.slx`| Simulink model of the analog reference circuit.        |
| `ssc_script.m`  | Matlab script that computes the sweeptone input signal.
| `ssc_output.mat`       | Ground-truth output from Simscape simulation.         |
| `Report.pdf`   | The assignment [report](Report.pdf) with detailed explanations and results. |

## Technical Implementation Details

### 1. WDF Architecture

- **Root Element:**   
The root of the WDF is the only non-adaptable element: the ideal voltage source 𝐹𝑖𝑛 (corresponding to a force generator in the mechanical domain in impedance analogy).

- **Series Junctions:**  
Denoted as **S1**, **S2** (the only five-port series junction; its scattering matrix was obtained using the fundamental loop matrix **B** method in order to optimize computational efficiency and
improve performances), **S3** and **S4** (connected to the resistor **R3**, from which the output
signal is extracted for comparison).

- **Parallel Junctions:**  
Denoted as **P1**, **P2**, **P3**

### 2. Configuration of Free Parameters 

The **adaptation conditions** for the individual ports are defined as:

$$ 
Z_2 = R_1, \quad Z_6 = \frac{T_s}{2C_1}, \quad Z_9 = R_2, \quad Z_{10} = \frac{2L_1}{T_s}, \quad Z_{11} = \frac{T_s}{2C_2},
$$
$$
Z_{14} = \frac{T_s}{2C_3},\quad Z_{17} = \frac{2L_2}{T_s}, \quad Z_{20} = \frac{T_s}{2C_4}, \quad Z_{22} = R_3, \quad Z_{23} = \frac{2L_3}{T_s}.
$$

In order to achieve **reflection-free adaptors**, the following composite impedances were computed:

$$
Z_{21} = Z_{22} + Z_{23}, \quad Z_{18} = \frac{Z_{19} \cdot Z_{20}}{Z_{19} + Z_{20}}, \quad Z_{15} = Z_{16} + Z_{17},
$$
$$ Z_{12} = \frac{Z_{13} \cdot Z_{14}}{Z_{13} + Z_{14}},\quad Z_7 = \sum_{i=8}^{11} Z_i, \quad Z_4 = \frac{Z_5 \cdot Z_6}{Z_5 + Z_6}, \quad Z_3 = Z_1 + Z_2.
$$

The expressions of the Scattering Matrices can be found in the code implementation.

## Results and Conclusions

For a detailed analysis of the final results, including a comparison between the ground truth and the simulation outcomes, please refer to the [report](Report.pdf).

This project was developed through a collaborative effort by [Matteo Di Giovanni](https://github.com/matteodigii) and [Alessandro Mancuso](https://github.com/AleMancusoPOLI).


