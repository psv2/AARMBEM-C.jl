# Two identical mirrored collinear 250 atom-long carbyne wires (vacuum)

This page links to examples of input files needed to compute thermal
energy transport powers in a system of compact molecular bodies. The
particular system under consideration is of two identical mirrored
collinear 250 atom-long carbyne wires (carbon chains) in vacuum at a
separation ``d`` from each other.

## Molecular configuration files

The input data files are as follows.

``B_{\mathrm{e}}`` for carbyne: [C250carbyneBe.txt](C250carbyneBe.txt)

``B_{\mathrm{I}}`` for carbyne: [C250carbyneBI.txt](C250carbyneBI.txt)

``K_{\mathrm{e}}`` for carbyne: [C250carbyneKe.txt](C250carbyneKe.txt)

``K^{(0)}_{\mathrm{I}}`` for carbyne: [C250carbyneKI_mod.txt](C250carbyneKI_mod.txt)

``M_{\mathrm{e}}`` for carbyne: [C250carbyneMe.txt](C250carbyneMe.txt)

``M_{\mathrm{I}}`` for carbyne: [C250carbyneMI.txt](C250carbyneMI.txt)

``Q_{\mathrm{e}}`` for carbyne: [C250carbyneQe.txt](C250carbyneQe.txt)

Atomic coordinates ``\{\vec{r}_{p}\}`` for carbyne: [C250carbynexyz.txt](C250carbynexyz.txt)

Matrix for ``k_{0}`` values for two collinear carbyne wires: [twoC250carbynecollinearDeltaKI_k0.txt](twoC250carbynecollinearDeltaKI_k0.txt)

Matrix for ``a_{0}`` values for two collinear carbyne wires: [twoC250carbynecollinearDeltaKI_a0.txt](twoC250carbynecollinearDeltaKI_a0.txt)

Matrix for ``E_{\infty}`` values for two collinear carbyne wires: [twoC250carbynecollinearDeltaKI_Einf.txt](twoC250carbynecollinearDeltaKI_Einf.txt)

Overall configuration file: [twoC250carbynecollinearconfig.txt](twoC250carbynecollinearconfig.txt)

## Transformation file

Transformation file: [twoC250carbynecollineartrans.txt](twoC250carbynecollineartrans.txt)

## Frequency file

Frequency file: [freq31.txt](freq31.txt)

(Many more frequencies must be used to replicate the output for vdW
interaction free energies and thermal radiation powers: these can in
principle be extracted from the corresponding output files themselves,
but this is hampered in practice by the truncation of the output files
in this example for space reasons. For clarity, in the example output
files used, for thermal energy transport, the input frequencies
``\omega = |w|`` used were ``w_{n} = n \times 10^{12}~\mathrm{rad/s}``
for ``n \in \{1, 2, \ldots, 500\}``.)

## Output files

Thermal energy transport powers with only radiative couplings: [radiation_twoC250carbynecollinearVACheat.out](radiation_twoC250carbynecollinearVACheat.out)

Thermal energy transport powers with only conductive couplings: [conduction_twoC250carbynecollinearVACheat.out](conduction_twoC250carbynecollinearVACheat.out)

Thermal energy transport powers with both couplings: [both_twoC250carbynecollinearVACheat.out](both_twoC250carbynecollinearVACheat.out)

!!! note

    The ordering of output lines may vary from one run to another.

## Standard commands to yield these output files

Thermal energy transport powers:
```
julia -p 1 AARMBEMC_heat.jl mollistfilename=twoC250carbynecollinearconfig.txt outfilename=twoC250carbynecollinearVACheat.out translistfilename=twoC250carbynecollineartrans.txt freqlistfilename=freq31.txt Genv=VAC
```

## Other possible examples for commands (not exhaustive)

Thermal radiation powers in the absence of EM retardation:
```
julia -p 1 AARMBEMC_heat.jl mollistfilename=twoC250carbynecollinearconfig.txt outfilename=twoC250carbynecollinearVACheat.out translistfilename=twoC250carbynecollineartrans.txt freqlistfilename=freq31.txt Genv=VAC nonretarded
```
The output file names will actually be
`radiation_nonretarded_twoC250carbynecollinearVACheat.out`,
`conduction_nonretarded_twoC250carbynecollinearVACheat.out`, and
`both_nonretarded_twoC250carbynecollinearVACheat.out`.
