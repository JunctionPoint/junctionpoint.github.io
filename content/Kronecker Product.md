---
title: Kronecker Product
date: 2023-07-29
---
For two matrices $A \in R^{n \times p}$ and $B \in R^{m \times q}$, the $mn \times pq$ matrix

$$
A \otimes B = 
\begin{bmatrix}
	a_{11}B & a_{12}B & \cdots & a_{1p}B \\
	a_{21}B & a_{22}B & \cdots & a_{2p}B \\
	\vdots & \vdots & \ddots & \vdots \\
	a_{n1}B & a_{n2}B & \cdots & a_{np}B \\
\end{bmatrix}
$$

is called the Kronecker product of $A$ and $B$[^kronecker].

[^kronecker]: Todd K Moon and Wynn C Stirling. Mathematical methods and algorithms for signal processing, volume 1, chapter The Kronecker product and Kroneeker sum. Prentice hall Upper Saddle River, NJ, 2000.
