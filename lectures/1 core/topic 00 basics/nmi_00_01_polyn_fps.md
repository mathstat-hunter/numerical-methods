# 00.01 basics: polynomials & floating point

<details>

<summary>code, setup and dependencues</summary>

```python
if True: # settings for the people
  import matplotlib.pyplot as plt
  import numpy as np
```

```python
def img_fps():
  es = [-2,-1,0,1,2]
  ms = np.arange(1,2,0.125) # implicit 1 + 3 position binary fractional
  zeros = np.zeros_like(ms)

  # plot
  plt.close("all")
  plt.figure(figsize=(20,2))

  plt.plot([0,8],[0,0],color="0.5")
  for e in es:
    xs = ms*pow(2,e)  # + mantissa x base^e
    plt.scatter(xs,zeros,marker="|",s=500)
    s_e = "$2^{" + f"{e}" + "}$"
    plt.text(xs.min(),.015,s_e,size="xx-large")

  plt.axis("off")
  plt.show()
```

</details>

## 0 intro

<b>numerical methods</b>, as distinguished from other branches of mathematics and from computer science,

1. work with arbitrary real numbers (including rational <b>approximations</b> of irrational numbers) and
2. consider <b>cost</b> and
3. consider <b>accuracy</b>.[^1]

<i>this class will provide another way to express, to extend your math.</i>

numerical methods are the algorithms; <b>numerical analysis</b> is the study of their properties -- ie, accuracy, stability, convergence, efficiency, usw.

## 1 polynomials

<i>The most fundamental operations of arithmetic are <b>addition</b> and <b>multiplication</b>. These are also the operations needed to evaluate a polynomial $p(x)$ at a particular value $x$. It is no coincidence that polynomials are the basic building blocks for many computational techniques we will construct.[^2]</i>

### i) evaluation

<details>

<summary>example 01</summary><br/>

consider $\enspace p(x) = a_4x^4 + a_3x^3 + a_2x^2 + a_1x + a_0$.

with computational considerations:

1. <b>approximate</b> $p(x)$ at $x$ while
2. minimizing <b>operations</b> and
3. maximizing <b>accuracy</b>.

wrt operations,

- method 1, step individually:
  - $p(x) = a_4 \times x \times x \times x \times x + a_3 \times x \times x \times x + a_2 \times x \times x + a_1 \times x + a_0 \mapsto 14$ operations.

- method 2, cache and reuse:
  - $x_2 = x \times x, x_3 = x_2 \times x, x_4 = x_3 \times x \mapsto 3$ operations;
  - $p_4 = a_4 \times x_4, p_3 = a_3 \times x_3, p_2 = a_2 \times x_2, p_1 = a_1 \times x_1 \mapsto 4$ operations;
  - $p(x) = p_4 + p_3 + p_2 + p_1 + a_0 \mapsto 4$ operations $\mapsto 11$ operations total.

- method 3, nested multiplication ([horners method](https://en.wikipedia.org/wiki/Horner%27s_method)):
  - $p(x) = (((a_4 \times x + a_3) \times x + a_2) \times x + a_1) \times x + a_0 \mapsto 8$ operations.

</details>

## 2 binary notation

### i) conversion to decimal

$$
\Rightarrow \enspace \dots b_2 \times 2^2 + b_1 \times 2^1 + b_0 \times 2^0 + b_{-1} \times 2^{-1} + b_{-2} \times 2^{-2} \dots
$$

<details>

<summary>example 02</summary><br/>

evaluate $111.11_2$.

$$
\begin{align}
  \text{integer:} &\quad 1 \times 2^2 + 1 \times 2^1 + 1 \times 2^0 = 4 + 2 + 1 = 7 \\
  \\
  \text{fractional:} &\quad 1 \times 2^{-1} + 1 \times 2^{-2} = \tfrac{1}{2} + \tfrac{1}{4} = \tfrac{3}{4} \\
  \\
  &\quad\Downarrow \\
  \\
  111.11_2 &= 7_{10} + (\tfrac{3}{4})_{10} = 7.75_{10}.
\end{align}
$$

</details>

### ii) conversion from decimal

<details>

<summary>example 03</summary><br/>

evaluate $111.25_{10}$.

$$
\begin{align}
  \text{integer:} &\qquad\; \tfrac{111}{2} = 55\, R\, 1 \\
  &\quad\rightarrow \tfrac{55}{2} = 27\, R\, 1 \\
  &\quad\rightarrow \tfrac{27}{2} = 13\, R\, 1 \\
  &\quad\rightarrow \tfrac{13}{2} = 6\, R\, 1 \\
  &\quad\rightarrow \;\tfrac{6}{2} = 3\, R\, 0 \\
  &\quad\rightarrow \;\tfrac{3}{2} = 1\, R\, 1 \\
  &\quad\rightarrow \;\tfrac{1}{2} = 0\, R\, 1 \\
  \\
  &\quad\rightarrow 1101111, \enspace\text{remainders in reverse order}  \\
  \\
  \text{fractional:} &\qquad\enspace 0.25\times 2 = 0.50 + 0 \\
  &\quad\rightarrow 0.50\times 2 = 0.00 + 1 \\
  \\
  &\quad\rightarrow 0.01, \enspace\text{integers in order from left to right} \\
  \\
  &\quad\Downarrow \\
  \\
  111.25_{10} &= 1101111_2 + 0.01_2 = 1101111.01_2.
\end{align}
$$

</details>

## 3 polynomials in the machine

### i) digital representation

$$
\begin{align}
  x &= [d_{N-1},\dots,d_1,d_0] \quad\text{digital vector} \\
  \\
  &= d_{N-1} \times b^{N-1} + \dots + d_1 \times b^1 + d_0\times b^0 \quad\text{with}\textbf{ precision } N \text{ and}\textbf{ base } b.
\end{align}
$$

<details>

<summary>example 04</summary><br/>

- base 10: $\quad 500_{10} = [5,0,0]; \quad [5] = 5_{10}$.
- base 02: $\quad [1,0,1] = 101_2 = 1\times 2^2 + 0\times 2^1 + 1\times 2^0 = 4 + 0 + 1 = 5_{10}$.

</details>

### ii) fixed/positional representation

<details>

<summary>example 04, continued</summary><br/>

- base 02: $101_\color{blue}{2} = 1\times \color{blue}{2}^2 + 0\times \color{blue}{2}^1 + 1\times \color{blue}{2}^0$

where RHS is <b>fixed representation</b> and LH subscript is the base or <b>radix</b> r.

</details>

additionally, precision $N\ge 1, r\ge 2$ such that

$$
\begin{align}
  x = &\enspace \sum^N d_kr^k \text{ has } r^N \textbf{ permutations}
\end{align}
$$

and can also be written as

$$
r^N = \color{blue}{(r-1)}({\bf r^{N-1}}) + \color{red}{(r^{N-1})} = \color{blue}{[r-1]_{N-1}}[r]_{N-2}\dots[r]_1[r]_0 + \color{red}{[r]_{-1}[r]_{-2}\dots [r]_{N-2}[r]_{N-1}}
$$

where subscripts denote position wrt exponent.

## resources

- horners method [@wiki](https://en.wikipedia.org/wiki/Horner%27s_method)
- telescoping sum [@wiki](https://en.wikipedia.org/wiki/Telescoping_series)
- floating-point [@wiki](https://en.wikipedia.org/wiki/Floating-point_arithmetic) [@youtube #1](https://www.youtube.com/watch?v=dQhj5RGtag0) [#2-pt1](https://www.youtube.com/watch?v=gc1Nl3mmCuY) [#2-pt2](https://www.youtube.com/watch?v=b2FgF2sUoS8)
- unit in last place (ulp) [@wiki](https://en.wikipedia.org/wiki/Unit_in_the_last_place)
- machine epsilon [@wiki](https://en.wikipedia.org/wiki/Machine_epsilon)
- [IEEE](https://www.ieee.org/) [754](https://en.wikipedia.org/wiki/IEEE_754)-[2019](https://standards.ieee.org/ieee/754/6210/)

## references

[^1]: johnson, sg. <i>[18.335, introduction to numerical methods](https://ocw.mit.edu/courses/18-335j-introduction-to-numerical-methods-spring-2019/),</i> mit.ocw, spring 2015.
[^2]: sauer, tim. <i>numerical analysis, 2nd edition</i>, pearson education, 2012, p1.
[^3]: martinez, vincent. <i>math 685</i>, hunter, spring 2023.
[^4]: <i>ibid</i>.
[^5]: nerdfirst. <i>[denormal numbers](https://www.youtube.com/watch?v=b2FgF2sUoS8)</i>, [0612 tv](https://www.youtube.com/@NERDfirst), 2020.
