// Return the truncation of Gaussian hypergeometric series G(a,b,c|z) by degree d.
function truncated_gauss(p,a,b,c,d)
    K := GF(p); T<z> := PolynomialRing(K);
    poch_a := 1; poch_b := 1; poch_c := 1; poch_1 := 1; sum := 1;
    for n in [1..d] do
        poch_a *:= a+n-1; poch_b *:= b+n-1; poch_c *:= c+n-1; poch_1 *:= n;
        sum +:= (poch_a*poch_b)/(poch_c*poch_1)*z^n;
    end for;
    return sum;
end function;

p := 13;
while p le 20000 do
    p := NextPrime(p); start := Realtime();
    if p mod 8 eq 1 then
        g1 := truncated_gauss(p,1/2,1/4,3/4,(p-1)/4);
        g2 := truncated_gauss(p,3/4,1/8,3/8,(p-1)/8);
        g3 := truncated_gauss(p,3/4,3/8,5/8,(3*p-3)/8);
        g4 := truncated_gauss(p,1/4,1/8,7/8,(p-1)/8);
        g5 := truncated_gauss(p,3/4,9/8,11/8,(p-9)/8);
        g6 := truncated_gauss(p,3/4,11/8,13/8,(3*p-11)/8);
        G := GCD([g1,g2,g3,g4,g5,g6]);
    elif p mod 8 eq 3 then
        g1 := truncated_gauss(p,1/2,3/4,5/4,(p-3)/4);
        g2 := truncated_gauss(p,3/4,3/8,5/8,(p-3)/8);
        g3 := truncated_gauss(p,1/4,1/8,7/8,(3*p-1)/8);
        g4 := truncated_gauss(p,1/4,3/8,9/8,(p-3)/8);
        g5 := truncated_gauss(p,3/4,11/8,13/8,(p-11)/8);
        G := GCD([g1,g2,g3,g4,g5]);
    elif p mod 8 eq 5 then
        g1 := truncated_gauss(p,1/2,1/4,3/4,(p-1)/4);
        g2 := truncated_gauss(p,1/4,-3/8,3/8,(p+3)/8);
        g3 := truncated_gauss(p,3/4,5/8,7/8,(p-5)/8);
        g4 := truncated_gauss(p,3/4,7/8,9/8,(3*p-7)/8);
        g5 := truncated_gauss(p,1/4,5/8,11/8,(p-5)/8);
        G := GCD([g1,g2,g3,g4,g5]);
    else
        g1 := truncated_gauss(p,1/2,3/4,5/4,(p-3)/4);
        g2 := truncated_gauss(p,3/4,7/8,9/8,(p-7)/8);
        G := GCD([g1,g2]);
    end if;
    if p mod 8 ne 7 then
        assert Degree(G) eq 0;  // expectation
    else
        printf "Characteristic = %5o:  Number = ", p;
        if p mod 16 eq 15 then (Degree(G)-1)/2; else Degree(G)/2; end if;
    end if;
end while;