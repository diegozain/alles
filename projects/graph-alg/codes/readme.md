# 🐊 error correcting codes

diego domenzain

2.2025

---

```text
   Hamming code (7,4)
 
    ⟶ correct one bit-flip for a word with 4 letters in 𝔽₂. 
 
   let P be the points of PG(2,2) that are not boring (🔴) as columns.
 
 
             🟣
             /|\
            / | \
           /  |  \
         🔴   |   🔴
         /    |    \
        /    🔴     \
       /      |      \
      /       |       \
     /        |        \
   🟣--------🔴---------🟣
 
 

  H = [ P | I₃ ]
      
      .---. 
  G = | I₄|
      |-P |
      .---.
 
 
 let w be the data, a vector of size 4 with entries in 𝔽₂,
 
     c = G·w ............. code-word
         ~~~  some noise  ~~~
         ~~~ happens here ~~~
     i = H⋅c ............. i is a vector representing a 
                           number in binary, 
                           where the error happened.
                           (it's really the column whose position
                           in H is the bit that got corrupted)
     j = ∑ i[k]·2ᵏ
  c[j] = c[j] + 1 ........ c gets corrected at position j.
```

[![](projects/graph-alg/pics/hamming74.png)](ipynb/lincodes.ipynb)