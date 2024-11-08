# EnerFlex Project

## Overview
EnerFlex is a comprehensive energy trading and management system designed to facilitate energy transactions and optimize operations within local energy communities. This project provides models and tools for peer-to-peer (P2P) energy trading, considering grid interactions and community-level energy management.

## Project Structure
```
EnerFlex/
|-- data/                    # Directory for datasets used in simulations
|-- models/                  # Core models and algorithms for energy trading
|-- scripts/                 # Script files for data processing and running simulations
|-- results/                 # Directory for storing simulation outputs
|-- README.md                # Project documentation
|-- requirements.txt         # Dependencies and libraries
|-- main.m                   # Main script for running the project (MATLAB)
```

## Features
- **Day-Ahead and Real-Time Market Operations**
- **P2P Energy Trading** between prosumers and consumers
- **Integration with Wholesale Market** for optimal energy trading
- **Battery Management System** for energy storage and utilization
- **Optimization Algorithms** using MATLAB and YALMIP

## Prerequisites
- MATLAB (R2020b or later)
- YALMIP toolbox
- CPLEX or Gurobi solver (optional for better performance)

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/rathorsumit2006/EnerFlex.git
   ```
2. Navigate to the project directory:
   ```bash
   cd EnerFlex
   ```
3. Install necessary dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage
1. **Initialize MATLAB Environment**:
   Ensure you have the required toolboxes installed.

2. **Run the Main Script**:
   ```matlab
   run('main.m')
   ```

3. **Configure Parameters**:
   Adjust input parameters in the `main.m` script or supporting configuration files to customize simulations.

## Configuration
- The project uses `.m` files for model implementation.
- Ensure the `data` directory contains all necessary input datasets.
- Modify `config.m` for different scenarios and parameter tuning.

## Contributing
We welcome contributions! Please follow the standard GitHub process:
1. Fork the repository.
2. Create a new branch (`feature/your-feature`).
3. Commit your changes.
4. Push to your branch and create a Pull Request.

## License
This project is licensed under the MIT License. See `LICENSE` for more information.

## Contact
For any questions or collaboration inquiries, please contact:
- **Rathor Sumitkumar K, Ph.D.**
- [rathorsumit2006@gmail.com](mailto:rathorsumit2006@gmail.com)
- [LinkedIn](https://www.linkedin.com/in/rathorsumit2010/)

---

Feel free to use or modify this project as needed, and star the repository if you find it helpful!
