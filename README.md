# Stabilizer Rényi Entropy

This repository provides an implementation to compute the [Stabilizer Rényi Entropy (SRE)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.128.050402) for qubit systems. Using a combination of linear algebra operations and Pauli string manipulations, the code calculates the \( n \)-th Rényi entropy for stabilizer states.

## Overview

The main computation involves evaluating expectation values of Pauli strings acting on spin configurations, then summing these contributions to obtain the Rényi entropy for a given quantum state. The code is parallel-ready.

The algorithm computes the \( n \)-th order Stabilizer Rényi Entropy (SRE) for a quantum state represented by a wavefunction \( \psi \) on \( N \) qubits. The SRE is a measure of entanglement specifically tailored for stabilizer states, and is defined as:

$$
S_n = \frac{1}{1 - n} \log \left( \sum_{P \in \mathcal{P}_N} \left| \langle \psi | P | \psi \rangle \right|^{2n} \right),
$$

we make use of the fact that 

$$
\langle \psi | P | \psi \rangle = \sum_s \overline{\psi(s')} \cdot \alpha \cdot \psi(s)
$$

to avoid ever building the matrices and having memory issues.

where:
- $S_n$ is the $n$-th order SRE,
- $\mathcal{P}_N$ is the set of all Pauli strings of length $N$,
- $P$ is a Pauli string, which is an operator of the form $P = P_1 \otimes P_2 \otimes \dots \otimes P_N$ with $P_i \in \{ \mathbb{I}, X, Y, Z \}$,
- $\langle \psi | P | \psi \rangle$ is the expectation value of the Pauli string $P$ on the state $|\psi\rangle$.

### Steps in the Algorithm:

1. **Initialize**: Set up a sum $\( \text{expectation} \textunderscore \text{sum} = 0 \)$.

2. **Iterate Over Pauli Strings**: For each Pauli string $P$ (i.e., each combination of Pauli operators across $N$ qubits),
   - Apply $P$ to the quantum state by computing the action of each Pauli operator on the individual qubit spins, obtaining a new spin configuration $s'$ and a complex coefficient $\alpha$.
   
3. **Compute Expectation Value for Each Pauli String**: For each spin configuration $s$,
   - Calculate: $\langle \psi | P | \psi \rangle = \sum_s \overline{\psi(s')} \cdot \alpha \cdot \psi(s)$, where $s'$ is the transformed spin configuration under $P$, and $\overline{\psi(s')}$ is the conjugate amplitude of the transformed configuration.

4. **Accumulate the \( 2n \)-th Power**: Compute \( \left| \langle \psi | P | \psi \rangle \right|^{2n} \) and add it to the `expectation_sum`.

5. **Final Computation of $S_n$**: After iterating over all Pauli strings,
   - Compute $S_n = \frac{1}{1 - n} \log \left( \frac{\text{expectation} \textunderscore \text{sum}}{2^N} \right)$, where $2^N$ is the normalization factor for the number of possible spin configurations.



### Usage

1. **Dependencies**: Make sure to have `LinearAlgebra`, `ProgressMeter`, `Distributed`, and `SharedArrays` installed.
2. **Running the Code**:
   - Define a quantum state vector `psi` as input.
   - Call `Mn(n, psi)` with your desired value of \( n \) and state vector `psi`.
   - The function will output the Stabilizer Rényi Entropy.

### Example

```julia
using LinearAlgebra, ProgressMeter

# Define a state vector `psi`
psi = [1/sqrt(2), 0, 0, 1/sqrt(2)]  # Example state (|00⟩ + |11⟩)/√2

# Compute the 2nd order Stabilizer Rényi Entropy
entropy = Mn(2, psi)
println("Stabilizer Rényi Entropy (n=2): ", entropy)

```

## Disclaimer

This code was collaboratively developed with @laurinbrunner. While we've focused on creating a functional implementation, there may be opportunities to optimize its performance further. We noticed a lack of similar implementations available, so we hope this repository helps fill a gap in the community. 

Contributions and suggestions are very welcome! If you have ideas for improvements or new features, feel free to submit a pull request.



