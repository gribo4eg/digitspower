using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;

namespace MultInverse
{
    static class Primality_Tests
    {  
     
        public static bool Prime_Test_Miller_Rabin(BigInteger y)
        {
            BigInteger b = 0, T, a;
            bool flag = true;
            int k = 0, i=0; //a - osnova
            if ((y & 1) == 0) return false;  // nuzhno ne4etnoe 4islo
            if (y == 3) return true;
            do
            {
                int N = (int)Math.Ceiling(Math.Ceiling(BigInteger.Log(y - 1, 2)) / 8); //kol-vo byte v 4isle (y-1)
                int rnd = Functions.rand(1, N);                //
                a = Functions.random_max(rnd);          //generiruem osnovu <<a>> v diapazone ot 1 do (y-1)
            }
            while ((a >= y - 1) || (a <= 1));

            
                if (BigInteger.GreatestCommonDivisor(a, y) != 1) return false;
                else
                {
                    while (flag == true)
                    {
                        k++;
                        b = (y - 1) >> k;                    // y-1 = 2^k * b,
                        if ((b & 1) == 1) flag = false;      // b - ne4etnoe
                    }
                    T = BigInteger.ModPow(a, b, y);
                    if ((T == 1) || (T == y - 1)) return true;
                       
                    else
                    {
                        for (i = 1; i < k; i++)
                        {
                            T = BigInteger.ModPow(T, 2, y);
                            if (T == y - 1) return true;
                        }
                    }
                }              
            
           return false;
        }
    }
}
