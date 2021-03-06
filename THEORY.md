
Theory
======

This part is still a draft. Please, see the McQuarrie "Statistical Mechanics" 
book [@mcquarrie2000statistical] for more information.

Partition Function
------------------


\begin{align}
Q &= \frac{q^N}{N!} \\
\ln Q &= N\ln q - N\ln N +N = N\ln \frac{qe}{N} \\
q &= q_{tr}\;q_{rot}\;q_{vib}\;q_{elec}
\end{align}

Translational degrees of freedom:
\begin{align}
q_{tr} &= \left( \frac{2 \pi m k_B T}{h^2}\right)^{t/2} V \\
\ln q_{tr} &= \frac{t}{2} \ln \left( \frac{2 \pi m k_B T}{h^2}\right) + \ln V \\
\frac{\partial}{\partial T} \ln q_{tr} &= \frac{t}{2} \frac{1}{T}
\end{align}

Rotational:
\begin{align}
q_{rot} &= \frac{\sqrt{\pi}}{\sigma} \left( \frac{8 \pi^2 k_B T }{h^2} \right)^{r/2} \left( \prod_{i=1}^r I_i \right)^{1/2} \\
\ln q_{rot} &= \frac{1}{2} \ln \frac{\pi}{\sigma^2} + \frac{r}{2}\ln\frac{8\pi^2k_BT}{h^2} + \frac{1}{2}\sum_{i=1}^r \ln I_i \\
\frac{\partial}{\partial T} \ln q_{rot} &= \frac{r}{2} \frac{1}{T}
\end{align}

Vibrational (classical limit):
\begin{align}
q_{vib,cl} &= \prod _{i=1} ^f \frac{k_BT}{h\nu_i} \\
\ln q_{vib,cl} &= \sum _{i=1} ^f \ln \frac{k_BT}{h\nu_i} \\
\frac{\partial}{\partial T} \ln q_{vib,cl} &= f \frac{1}{T}
\end{align}

Vibrational (quantum):
\begin{align}
q_{vib,qm} &= \prod _{i=1} ^f \frac{\exp\frac{h\nu}{2k_BT}}{1-\exp\frac{h\nu}{k_BT}} = \prod _{i=1} ^f \left[ 2\sinh \left(\frac{h\nu}{2k_BT}\right)\right] ^{-1} \\
\ln q_{vib,qm} &= -\sum _{i=1} ^f \ln \left[ 2\sinh \left(\frac{h\nu}{2k_BT}\right)\right] \\
\frac{\partial}{\partial T} \ln q_{vib,qm} &= \frac{1}{T} \sum _{i=1} ^f \frac{\frac{h\nu}{2k_BT}}{\tanh\left(\frac{h\nu}{2k_BT}\right)}
\end{align}

Electronic:
\begin{align}
q_{elec} &= \exp \left(-\frac{E_m}{k_BT}\right) \\
\ln q_{elec} &= -\frac{E_m}{k_BT} \\
\frac{\partial}{\partial T} \ln q_{elec} &= \frac{E_m}{k_BT} \frac{1}{T}
\end{align}


Helmholtz Free energy
----------------------

\begin{align}
F = U-TS
\end{align}

\begin{align}
F &= -k_BT\ln Q = -Nk_BT\ln\frac{qe}{N} \\
  &= -Nk_BT\ln\frac{q_{tr}e}{N} - Nk_BT\ln q_{rot} - Nk_BT\ln q_{vib} - Nk_BT\ln q_{elec}\\
  &= F_{tr} + F_{rot} + F_{vib} + F_{elec}
\end{align}

\begin{align}
F_{tr} &= -Nk_BT\ln\frac{q_{tr}e}{N} = -Nk_BT\ln\left[ \left( \frac{2 \pi m k_B T}{h^2}\right)^{t/2} \frac{Ve}{N}\right] \\
       &= -Nk_BT \left( \frac{t}{2}\ln \frac{2 \pi m k_B T}{h^2} + 1 - \ln\frac{N}{V} \right)
\end{align}
\begin{align}
F_{rot} = -Nk_BT\ln q_{rot} = -Nk_BT \left[ \frac{1}{2} \ln \frac{\pi}{\sigma^2} + \frac{r}{2}\ln\frac{8\pi^2k_BT}{h^2} + \frac{1}{2}\sum_{i=1}^r \ln I_i\right]
\end{align}
\begin{align}
F_{vib,cl} = -Nk_BT\ln q_{vib,cl} = -Nk_BT \sum_{i=1}^f \ln \frac{k_BT}{h\nu_i}
\end{align}
\begin{align}
F_{vib,qm} = -Nk_BT\ln q_{vib,qm} = Nk_BT \sum _{i=1} ^f \ln \left[ 2\sinh \left(\frac{h\nu}{2k_BT}\right)\right]
\end{align}
\begin{align}
F_{elec} = -Nk_BT\ln q_{elec} = NE_m
\end{align}


Molar Helmholtz free energy
---------------------------

\begin{align}
\mu &= F_m = \frac{\partial F}{\partial N} \\
    &= \frac{\partial F_{tr}}{\partial N} + \frac{\partial F_{rot}}{\partial N} + \frac{\partial F_{vib}}{\partial N} + \frac{\partial F_{elec}}{\partial N} \\
    &= \mu_{tr} + \mu_{rot} + \mu_{vib} + \mu_{elec}
\end{align}
\begin{align}
\mu_{tr} &= \frac{\partial F_{tr}}{\partial N} = -k_BT\frac{t}{2}\ln \frac{2 \pi m k_B T}{h^2} + k_BT\ln\frac{N}{V} = \frac{F_{tr}}{N} + k_BT
\end{align}
\begin{align}
\mu_{rot} &= \frac{\partial F_{rot}}{\partial N} = -k_BT \left[ \frac{1}{2} \ln \frac{\pi}{\sigma^2} + \frac{r}{2}\ln\frac{8\pi^2k_BT}{h^2} + \frac{1}{2}\sum_{i=1}^r \ln I_i\right] = \frac{F_{rot}}{N}
\end{align}
\begin{align}
\mu_{vib,cl} &= \frac{\partial F_{vib,cl}}{\partial N} = -k_BT \sum_{i=1}^f \ln \frac{k_BT}{h\nu_i} = \frac{F_{vib,cl}}{N}
\end{align}
\begin{align}
\mu_{vib,qm} &= \frac{\partial F_{vib,qm}}{\partial N} = k_BT \sum _{i=1} ^f \ln \left[ 2\sinh \left(\frac{h\nu}{2k_BT}\right)\right] = \frac{F_{vib,qm}}{N}
\end{align}
\begin{align}
\mu_{elec} &= \frac{\partial F_{elec}}{\partial N} = E_m = \frac{F_{elec}}{N}
\end{align}


Internal energy
------------------------
\begin{align}
U &= k_BT^2\frac{\partial}{\partial T}\ln Q = Nk_BT^2\frac{\partial}{\partial T}\ln q \\
  &= Nk_BT^2\frac{\partial}{\partial T}\ln q_{tr} + Nk_BT^2\frac{\partial}{\partial T}\ln q_{rot} + Nk_BT^2\frac{\partial}{\partial T}\ln q_{vib} + Nk_BT^2\frac{\partial}{\partial T}\ln q_{elec}\\
  &= U_{tr} + U_{rot} + U_{vib} + U_{elec}
\end{align}
\begin{align}
U_{tr} = Nk_BT^2\frac{\partial}{\partial T}\ln q_{tr} = Nk_BT\frac{t}{2}
\end{align}
\begin{align}
U_{rot} = Nk_BT^2\frac{\partial}{\partial T}\ln q_{rot} = Nk_BT\frac{r}{2}
\end{align}
\begin{align}
U_{vib,cl} = Nk_BT^2\frac{\partial}{\partial T}\ln q_{vib,cl} = Nk_BTf
\end{align}
\begin{align}
U_{vib,qm} = Nk_BT^2\frac{\partial}{\partial T}\ln q_{vib,qm} = Nk_BT \sum _{i=1} ^f \frac{\frac{h\nu}{2k_BT}}{\tanh\left(\frac{h\nu}{2k_BT}\right)}
\end{align}
\begin{align}
U_{elec} = Nk_BT^2\frac{\partial}{\partial T}\ln q_{elec} = NE_m
\end{align}

Molar internal energy
---------------------

\begin{align}
U_m &= \frac{\partial U}{\partial N} = \frac{\partial U_{tr}}{\partial N} + \frac{\partial U_{rot}}{\partial N} + \frac{\partial U_{vib}}{\partial N} + \frac{\partial U_{elec}}{\partial N} = U_{m,tr} + U_{m,rot} + U_{m,vib} + U_{m,elec}
\end{align}
\begin{align}
U_{m,tr} = \frac{\partial U_{tr}}{\partial N} = k_BT\frac{t}{2} = \frac{U_{tr}}{N}
\end{align}
\begin{align}
U_{m,rot} = \frac{\partial U_{rot}}{\partial N} = k_BT\frac{r}{2} = \frac{U_{rot}}{N}
\end{align}
\begin{align}
U_{m,vib,cl} = \frac{\partial U_{vib,cl}}{\partial N} = k_BTf = \frac{U_{vib,cl}}{N}
\end{align}
\begin{align}
U_{m,vib,qm} = \frac{\partial U_{vib,qm}}{\partial N} = k_BT \sum _{i=1} ^f \frac{\frac{h\nu}{2k_BT}}{\tanh\left(\frac{h\nu}{2k_BT}\right)} = \frac{U_{vib,qm}}{N}
\end{align}
\begin{align}
U_{m,elec} = \frac{\partial U_{elec}}{\partial N} = E_m = \frac{U_{elec}}{N}
\end{align}


Entropy
----------------
\begin{align}
S &= \frac{\partial}{\partial T}\left( k_BT\ln Q\right) = k_B\ln Q + k_BT\frac{\partial}{\partial T}\ln Q = Nk_B\ln\frac{qe}{N} + Nk_BT\frac{\partial}{\partial T}\ln q \\
  &= \left( Nk_B\ln \frac{q_{tr}e}{N} + Nk_BT\frac{\partial}{\partial T} \ln q_{tr} \right) 
   + \left( Nk_B\ln q_{rot} + Nk_BT\frac{\partial}{\partial T} \ln q_{rot} \right)  \nonumber \\ 
   &\qquad\qquad+\left( Nk_B\ln q_{vib} + Nk_BT\frac{\partial}{\partial T} \ln q_{vib} \right) 
   + \left( Nk_B\ln q_{elec} + Nk_BT\frac{\partial}{\partial T} \ln q_{elec} \right) \\
  &= S_{tr} + S_{rot} + S_{vib} + S_{elec}
\end{align}
\begin{align}
S_{tr} = Nk_B\ln \frac{q_{tr}e}{N} + Nk_BT\frac{\partial}{\partial T} \ln q_{tr} = Nk_B\left(\frac{t}{2}\ln\frac{2\pi m k_BT}{h^2} + \frac{t+2}{2} - \ln\frac{N}{V}\right)
\end{align}
\begin{align}
S_{rot} = Nk_B\ln q_{rot} + Nk_BT\frac{\partial}{\partial T} \ln q_{rot} = Nk_B\left( \frac{1}{2} \ln \frac{\pi}{\sigma^2} + \frac{r}{2}\ln\frac{8\pi^2k_BT}{h^2} + \frac{1}{2}\sum_{i=1}^r \ln I_i + \frac{r}{2} \right)
\end{align}
\begin{align}
S_{vib,cl} = Nk_B\ln q_{vib,cl} + Nk_BT\frac{\partial}{\partial T} \ln q_{vib,cl} = Nk_B \left(f+ \sum _{i=1} ^f \ln \frac{k_BT}{h\nu_i} \right) 
\end{align}
\begin{align}
S_{vib,qm} = Nk_B\ln q_{vib,qm} + Nk_BT\frac{\partial}{\partial T} \ln q_{vib,qm} = Nk_B \sum _{i=1} ^f \left[ \frac{\frac{h\nu}{2k_BT}}{\tanh\frac{h\nu}{2k_BT}}- \ln\left( 2\sinh\frac{h\nu}{2k_BT}\right) \right]
\end{align}
\begin{align}
S_{elec} = Nk_B\ln q_{elec} + Nk_BT\frac{\partial}{\partial T} \ln q_{elec} = 0
\end{align}

Molar Entropy
-----------------------
\begin{align}
S_m &= \frac{\partial S}{\partial N} = \frac{\partial S_{tr}}{\partial N} + \frac{\partial S_{rot}}{\partial N} + \frac{\partial S_{vib}}{\partial N} + \frac{\partial S_{elec}}{\partial N} = S_{m,tr} + S_{m,rot} + S_{m,vib} + S_{m,elec}
\end{align}
\begin{align}
S_{m,tr} = \frac{\partial S_{tr}}{\partial N} = k_B\left(\frac{t}{2}\ln\frac{2\pi m k_BT}{h^2} + \frac{t}{2} - \ln\frac{N}{V}\right) = \frac{S_{tr}}{N} - k_B
\end{align}
\begin{align}
S_{m,rot} = \frac{\partial S_{rot}}{\partial N} = k_B\left( \frac{1}{2} \ln \frac{\pi}{\sigma^2} + \frac{r}{2}\ln\frac{8\pi^2k_BT}{h^2} + \frac{1}{2}\sum_{i=1}^r \ln I_i + \frac{r}{2} \right) = \frac{S_{rot}}{N}
\end{align}
\begin{align}
S_{m,vib,cl} = \frac{\partial S_{vib,cl}}{\partial N} = k_B \left(f+ \sum _{i=1} ^f \ln \frac{k_BT}{h\nu_i} \right) = \frac{S_{vib,cl}}{N}
\end{align}
\begin{align}
S_{m,vib,qm} = \frac{\partial S_{vib,qm}}{\partial N} = k_B \sum _{i=1} ^f \left[ \frac{\frac{h\nu}{2k_BT}}{\tanh\frac{h\nu}{2k_BT}}- \ln\left( 2\sinh\frac{h\nu}{2k_BT}\right) \right] = \frac{S_{vib,qm}}{N}
\end{align}
\begin{align}
S_{m,elec} = \frac{\partial S_{elec}}{\partial N} = 0
\end{align}

