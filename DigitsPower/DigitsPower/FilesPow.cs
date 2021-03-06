using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Windows.Forms;
using System.Numerics;

namespace DigitsPower
{
    abstract class BinaryPow
    {
        public Inverse InvMethod { get; set; }
        public Multiply MultMethod { get; set; }

        public string Found { get; set; }
        public string Degree { get; set; }
        public string Mod { get; set; }
        public string Choice { get; set; }

        public BinaryPow(string found, string degree, string mod, string by)
        {
            Found = found;
            Degree = degree;
            Mod = mod;
            Choice = by;
        }

        public abstract string Name();

        protected abstract void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod);

        protected void SendStatus()
        {
            (Application.OpenForms[0] as MainForm).ChangeStatus(Name());
        }

        protected int[] MakeDigits(string[] mas)
        {
            int[] mas_i = new int[mas.Length]; ;

            for (int i = 0; i < mas.Length; i++)
            {
                var s = mas[i].Split('\\');
                mas_i[i] = Int32.Parse(s[s.Length - 1].Split('.')[0]);
            }
            Array.Sort(mas_i);
            return mas_i;
        }

        protected virtual void Gen(string[] one, string[] two, string[] three, string path, string di, string One, string Two, string Three)
        {
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);


            for (int f_len = 0; f_len < one.Length; f_len++)
            {
                lock (this)
                {
                    using (FileStream fin = new FileStream(di + "\\" + One.Split('\\')[0] + "_" + one_i[f_len] + ".csv", FileMode.Create, FileAccess.Write))
                    {
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write(Two.Split('\\')[0] + "\\" + Three.Split('\\')[0] + ";");

                            for (int m_len = 0; m_len < three.Length; m_len++)
                            {
                                sw.Write(three_i[m_len] + ";");
                            }
                            sw.WriteLine();
                            for (int d_len = 0; d_len < two.Length; d_len++)
                            {
                                sw.Write(two_i[d_len] + ";");
                                for (int m_len = 0; m_len < three.Length; m_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                   path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                   path + "\\" + Three + "\\" + three_i[m_len] + ".txt") + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                    }
                }
            }
        }

        public virtual void Create_Result()
        {
            SendStatus();
            string path = Directory.GetCurrentDirectory();
            string[] founds = Directory.GetFiles(path + "\\Base\\" + Found);
            string[] degrees = Directory.GetFiles(path + "\\Exponent\\" + Degree);
            string[] mods = Directory.GetFiles(path + "\\Modulus\\" + Mod);

            var mont = AdditionalParameters.montFlag ? "Mont_" : "";
            DirectoryInfo di;
            switch (Choice)
            {
                case "Base":
                    di = Directory.CreateDirectory($"{path}\\Results\\{mont}{Name()}_{Mod.Split('#')[0]}_{Found.Split('#')[0]}+_{Degree.Split('#')[0]}#{DateTime.Now.ToLocalTime().ToString().Replace(':', '-')}");
                    Gen(founds, degrees, mods, path, di.FullName, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                    break;
                case "Exponent":
                    di = Directory.CreateDirectory(String.Format("{0}\\{1}\\{7}{2}_{3}_{4}_{5}#{6}",
                                                   path, "Results", Name(), Mod.Split('#')[0], Found.Split('#')[0], 
                                                   Degree.Split('#')[0], DateTime.Now.ToLocalTime().ToString().Replace(':', '-'), mont));
                    Gen(degrees, founds, mods, path, di.FullName, "Exponent\\" + Degree, "Base\\" + Found, "Modulus\\" + Mod);
                    break;
                case "Modulus":
                    di = Directory.CreateDirectory(String.Format("{0}\\{1}\\{7}{2}_{3}+_{4}_{5}#{6}",
                                                   path, "Results", Name(), Mod.Split('#')[0], Found.Split('#')[0], 
                                                   Degree.Split('#')[0], DateTime.Now.ToLocalTime().ToString().Replace(':', '-'), mont));
                    Gen(mods, founds, degrees, path, di.FullName, "Modulus\\" + Mod, "Base\\" + Found, "Exponent\\" + Degree);
                    break;
                default:
                    try
                    {
                        di = Directory.CreateDirectory(String.Format("{0}\\{1}\\{7}{2}_{3}_{4}+_{5}#{6}",
                                                       path, "Results", Name(), Mod.Split('#')[0], Found.Split('#')[0], 
                                                       Degree.Split('#')[0], DateTime.Now.ToLocalTime().ToString().Replace(':', '-'), mont));
                        Gen(founds, degrees, mods, path, di.FullName, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                    }
                    catch (Exception err)
                    {
                        MessageBox.Show(err.Message, "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                    }
                    break;
            }

            SendStatus();
        }

        protected double Watch(string Found_file, string Degree_file, string Mod_file)
        {
            Stopwatch stw = new Stopwatch();
            List<double> sums = new List<double>();
            double f_sum = 0;
            int f_count = 0;
            
            using (StreamReader srf = new StreamReader(Found_file))
            {
                string found;
                while ((found = srf.ReadLine()) != null)
                {
                    double d_sum = 0;
                    int d_count = 0;
                    using (StreamReader srd = new StreamReader(Degree_file))
                    {
                        string degree;
                        while ((degree = srd.ReadLine()) != null)
                        {
                            double m_sum = 0;
                            int m_count = 0;
                            using (StreamReader srm = new StreamReader(Mod_file))
                            {
                                string mod;
                                while ((mod = srm.ReadLine()) != null)
                                {
                                    stw.Start();
                                    LoopFunc(BigInteger.Parse(found), BigInteger.Parse(degree), BigInteger.Parse(mod));
                                    stw.Stop();

                                    m_sum +=  stw.Elapsed.TotalMilliseconds;
                                    m_count++;
                                    stw.Reset();
                                }

                            }
                            d_sum += m_sum / m_count;
                            d_count++;
                        }
                    }
                    f_sum += d_sum / d_count;
                    f_count++;
                }
            }
            return f_sum / f_count;
        }
    }

    abstract class BinaryPowUp
    {
        public string Found { get; set; }
        public string Degree { get; set; }
        public string Mod { get; set; }
        public string[] Choice { get; set; }
        List<int> aMax, bMax;

        public BinaryPowUp(string found, string degree, string mod, string by, string a_max, string b_max)
        {
            Found = found;
            Degree = degree;
            Mod = mod;
            Choice = by.Split(' ');
            aMax = GenFunctions.ReadString(a_max);
            bMax = GenFunctions.ReadString(b_max);
        }

        public abstract string Name();

        protected abstract void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, params BigInteger[] list);

        protected void SendStatus()
        {
            (Application.OpenForms[0] as MainForm).ChangeStatus(Name());
        }

        protected int[] MakeDigits(string[] mas)
        {
            int[] mas_i = new int[mas.Length];

            for (int i = 0; i < mas.Length; i++)
            {
                var s = mas[i].Split('\\');
                mas_i[i] = Int32.Parse(s[s.Length - 1].Split('.')[0]);
            }
            Array.Sort(mas_i);
            return mas_i;
        }

        private string AddZeros(string num, string maxNum)
        {
            while (num.Length < maxNum.Length)
            {
                num = num.Insert(0, "0");
            }
            return num;
        }

        private string AddZeros(int num, int maxNum)
        {
            return AddZeros(num.ToString(), maxNum.ToString());
        }

        protected List<BigInteger> GetFileDigits(string path)
        {
            List<BigInteger> result = new List<BigInteger>();
            using (StreamReader file = new StreamReader(path))
            {
                string digit;
                while ((digit = file.ReadLine()) != null)
                    result.Add(BigInteger.Parse(digit));
            }
            return result;
        }

        private DirectoryInfo CD(string[] by, string Found, string Degree, string Mod)
        {
            DirectoryInfo di;
            string path = Directory.GetCurrentDirectory();
            di = Directory.CreateDirectory(String.Format("{0}\\{1}\\{8}{2}_A({3}-{4})B({5}-{6})_(" + by[0][0] + "_" + by[1][0] + "_" + by[2][0] + ")#{7}",
                                         path, "Results", Name(), aMax[0], aMax[aMax.Count - 1],
                                         bMax[0], bMax[bMax.Count - 1],
                                         DateTime.Now.ToLocalTime().ToString().Replace(':', '-'), AdditionalParameters.montFlag ? "Mont_" : ""));

            using (StreamWriter sw = new StreamWriter(new FileStream(di.FullName + "\\" + "Info.txt", FileMode.Create, FileAccess.Write)))
            {
                sw.WriteLine(Mod.Split('#')[0]);
                sw.WriteLine(Found.Split('#')[0]);
                sw.WriteLine(Degree.Split('#')[0]);
            }

            return di;
        }

        #region GenMethods

        protected virtual void MaxBaseMaxGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three,  bool order = true)
         {
             string path = Directory.GetCurrentDirectory();
             int[] one_i = MakeDigits(one);
             int[] two_i = MakeDigits(two);
             int[] three_i = MakeDigits(three);
             DirectoryInfo amax_di, bmax_di;
             if (!order)
             {
                 var x = aMax;
                 aMax = bMax;
                 bMax = x;
             }
             for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
             {
                 amax_di = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                 for (int f_len = 0; f_len < one.Length; f_len++)
                 {
                     bmax_di = Directory.CreateDirectory(di + "\\" + amax_di.Name + "\\" + Choice[1] + "_" + AddZeros(one_i[f_len], one_i.Max()));
                     for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                     {
                         FileStream fin = new FileStream(di + "\\" + amax_di.Name + "\\" + bmax_di.Name + "\\" + Choice[2] + "_" + AddZeros(bMax[bmax_len], bMax.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                       using (StreamWriter sw = new StreamWriter(fin))
                       {
                           sw.Write(Two.Split('\\')[0] + "\\" + Three.Split('\\')[0] + ";");
 
                           for (int m_len = 0; m_len < three.Length; m_len++)
                           {
                               sw.Write(three_i[m_len] + ";");
                           }
                           sw.WriteLine();
                           for (int d_len = 0; d_len < two.Length; d_len++)
                           {
                               sw.Write(two_i[d_len] + ";");
                               for (int m_len = 0; m_len < three.Length; m_len++)
                               {
                                   sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                order ? aMax[amax_len] : bMax[bmax_len], 
                                                order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                               }
                               sw.WriteLine();
                           }
                       }
                       fin.Close();
                     }
                 }
             }
         }
 
         protected virtual void MaxExpMaxGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
         {
             string path = Directory.GetCurrentDirectory();
             int[] one_i = MakeDigits(one);
             int[] two_i = MakeDigits(two);
             int[] three_i = MakeDigits(three);
             DirectoryInfo amax_di, bmax_di;
             if (!order)
             {
                 var x = aMax;
                 aMax = bMax;
                 bMax = x;
             }
             for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
             {
                 amax_di = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                 for (int d_len = 0; d_len < two.Length; d_len++)
                 {
                     bmax_di = Directory.CreateDirectory(di + "\\" + amax_di.Name + "\\" + Choice[1] + "_" + AddZeros(two_i[d_len], two_i.Max()));
                     for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                     {
                         FileStream fin = new FileStream(di + "\\" + amax_di.Name + "\\" + bmax_di.Name + "\\" + Choice[2] + "_" + AddZeros(bMax[bmax_len], bMax.Max())  + ".csv", FileMode.Create, FileAccess.Write);
 
                         using (StreamWriter sw = new StreamWriter(fin))
                         {
                             sw.Write(One.Split('\\')[0] + "\\" + Three.Split('\\')[0] + ";");
 
                             for (int m_len = 0; m_len < three.Length; m_len++)
                             {
                                 sw.Write(three_i[m_len] + ";");
                             }
                             sw.WriteLine();
                             for (int f_len = 0; f_len < one.Length; f_len++)
                             {
                                 sw.Write(one_i[f_len] + ";");
                                 for (int m_len = 0; m_len < three.Length; m_len++)
                                 {
                                     sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                  path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                  path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                  order ? aMax[amax_len] : bMax[bmax_len], 
                                                  order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                 }
                                 sw.WriteLine();
                             }
                         }
                         fin.Close();
                     }
                 }
             }
         }
 
         protected virtual void MaxBaseExpGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
         {
             string path = Directory.GetCurrentDirectory();
             int[] one_i = MakeDigits(one);
             int[] two_i = MakeDigits(two);
             int[] three_i = MakeDigits(three);
             DirectoryInfo amax_di, bmax_di;
             if (!order)
             {
                 var x = aMax;
                 aMax = bMax;
                 bMax = x;
             }
             for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
             {
                 amax_di = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                 for (int f_len = 0; f_len < one.Length; f_len++)
                 {
                     bmax_di = Directory.CreateDirectory(di + "\\" + amax_di.Name + "\\" + Choice[1] + "_" + AddZeros(one_i[f_len], one_i.Max()));
                     for (int d_len = 0; d_len < two.Length; d_len++)
                     {
                         FileStream fin = new FileStream(di + "\\" + amax_di.Name + "\\" + bmax_di.Name + "\\" + Choice[2] + "_" + AddZeros(two_i[d_len], two_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                         using (StreamWriter sw = new StreamWriter(fin))
                         {
                             sw.Write(Three.Split('\\')[0] + "\\" + (order ? "b_max" : "a_max") + ";");
 
                             for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++) 
                             {
                                 sw.Write(bMax[bmax_len] + ";");
                             }
                             sw.WriteLine();
                             for (int m_len = 0; m_len < three.Length; m_len++)
                             {
                                 sw.Write(three_i[m_len] + ";");
                                 for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++) 
                                 {
                                     sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                  path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                  path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                  order ? aMax[amax_len] : bMax[bmax_len], 
                                                  order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                 }
                                 sw.WriteLine();
                             }
                         }
                         fin.Close();
                     }
                 }
             }
         }
        
        protected virtual void BaseMaxMaxGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo amax_di, one_dir;
            if (!order)
            {
                var x = aMax;
                aMax = bMax;
                bMax = x; 
            }
            for (int f_len = 0; f_len < one.Length; f_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(one_i[f_len], one_i.Max()));
                for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                {
                    amax_di = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                    for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + amax_di.Name + "\\" + Choice[2] + "_" + AddZeros(bMax[bmax_len], bMax.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write(Two.Split('\\')[0] + "\\" + Three.Split('\\')[0] + ";");
                            
                            for (int m_len = 0; m_len < three.Length; m_len++)
                            {
                                sw.Write(three_i[m_len] + ";");
                            }
                            sw.WriteLine();
                            for (int d_len = 0; d_len < two.Length; d_len++)
                            {
                                sw.Write(two_i[d_len] + ";");
                                for (int m_len = 0; m_len < three.Length; m_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                 path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                 path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                 order ? aMax[amax_len] : bMax[bmax_len], 
                                                 order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
 
         protected virtual void BaseExpMaxGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
         {
             string path = Directory.GetCurrentDirectory();
             int[] one_i = MakeDigits(one);
             int[] two_i = MakeDigits(two);
             int[] three_i = MakeDigits(three);
             DirectoryInfo one_dir, two_dir;
             if (!order)
             {
                 var x = aMax;
                 aMax = bMax;
                 bMax = x;
             }
             for (int f_len = 0; f_len < one.Length; f_len++)
             {
                 one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(one_i[f_len], one_i.Max()));
                 for (int d_len = 0; d_len < two.Length; d_len++)
                 {
                     two_dir = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(two_i[d_len], two_i.Max()));
                     for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                     {
                         FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + two_dir.Name + "\\" + Choice[2] + "_" + AddZeros(aMax[amax_len], aMax.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                         using (StreamWriter sw = new StreamWriter(fin))
                         {
                             sw.Write((order ? "b_max" : "a_max") + "\\" + Three.Split('\\')[0] + ";");
 
                             for (int m_len = 0; m_len < three.Length; m_len++)
                             {
                                 sw.Write(three_i[m_len] + ";");
                             }
                             sw.WriteLine();
                             for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                             {
                                 sw.Write(bMax[bmax_len] + ";");
                                 for (int m_len = 0; m_len < three.Length; m_len++)
                                 {
                                     sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                  path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                  path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                  order ? aMax[amax_len] : bMax[bmax_len], 
                                                  order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                 }
                                 sw.WriteLine();
                             }
                         }
                         fin.Close();
                     }
                 }
             }
         }
 
         protected virtual void BaseExpModGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
         {
             string path = Directory.GetCurrentDirectory();
             DirectoryInfo one_dir, two_dir;
             if (!order)
             {
                 var x = two;
                 two = three;
                 three = x;

                 var y = Two;
                 Two = Three;
                 Three = y;
             }
             int[] one_i = MakeDigits(one);
             int[] two_i = MakeDigits(two);
             int[] three_i = MakeDigits(three);
             for (int f_len = 0; f_len < one.Length; f_len++)
             {
                 one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(one_i[f_len], one_i.Max()));
                 for (int d_len = 0; d_len < two.Length; d_len++)
                 {
                     two_dir = Directory.CreateDirectory(di + "\\" + one_dir.Name +  "\\" + Choice[1] + "_" + AddZeros(two_i[d_len], two_i.Max()));
                     for  (int m_len = 0; m_len < three.Length; m_len++)
                     {
                         FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + two_dir.Name + "\\" + Choice[2] + "_" + AddZeros(three_i[m_len], three_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                         using (StreamWriter sw = new StreamWriter(fin))
                         {
                             sw.Write("b_max" + "\\" + "a_max;");
 
                             for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                             {
                                 sw.Write(aMax[amax_len] + ";");
                             }
                             sw.WriteLine();
                             for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                             {
                                 sw.Write(bMax[bmax_len] + ";");
                                 for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                                 {
                                     sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                  path + "\\" + 
                                                  (order ? Two + "\\" + two_i[d_len] : Three + "\\" + three_i[m_len]) + ".txt",
                                                  path + "\\" + 
                                                  (order ? Three + "\\" + three_i[m_len] : Two + "\\" + two_i[d_len]) + ".txt", 
                                                  aMax[amax_len], bMax[bmax_len]) + ";");
                                 }
                                 sw.WriteLine();
                             }
                         }
                         fin.Close();
                     }
                 }
             }
         }
        
        protected virtual void BaseModMaxGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
         {
             string path = Directory.GetCurrentDirectory();
             int[] one_i = MakeDigits(one);
             int[] two_i = MakeDigits(two);
             int[] three_i = MakeDigits(three);
             DirectoryInfo one_dir, two_dir;
             if (!order)
             {
                 var x = aMax;
                 aMax = bMax;
                 bMax = x;
             }
             for (int f_len = 0; f_len < one.Length; f_len++)
             {
                 one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(one_i[f_len], one_i.Max()));
                 for (int m_len = 0; m_len < three.Length; m_len++)
                 {
                     two_dir = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(three_i[m_len], three_i.Max()));
                     for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                     {
                         FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + two_dir.Name + "\\" + Choice[2] + "_" + AddZeros(aMax[amax_len], aMax.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                         using (StreamWriter sw = new StreamWriter(fin))
                         {
                             sw.Write((order ? "b_max" : "a_max") + "\\" + Two.Split('\\')[0] + ";");
 
                             for (int d_len = 0; d_len < two.Length; d_len++)
                             {
                                 sw.Write(two_i[d_len] + ";");
                             }
                             sw.WriteLine();
                             for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                             {
                                 sw.Write(bMax[bmax_len] + ";");
                                 for (int d_len = 0; d_len < two.Length; d_len++)
                                 {
                                     sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                  path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                  path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                  order ? aMax[amax_len] : bMax[bmax_len], 
                                                  order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                 }
                                 sw.WriteLine();
                             }
                         }
                         fin.Close();
                     }
                 }
             }
         }
        
        protected virtual void BaseMaxExpGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo amax_di, one_dir;
            if (!order)
            {
               var x = aMax;
                aMax = bMax;
                bMax = x;
            }
            for (int f_len = 0; f_len < one.Length; f_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(one_i[f_len], one_i.Max()));
                for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                {
                    amax_di = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                    for (int d_len = 0; d_len < two.Length; d_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + amax_di.Name + "\\" + Choice[2] + "_" + AddZeros(two_i[d_len], two_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write((order ? "b_max" : "a_max") + "\\" + Three.Split('\\')[0] + ";");
                            
                            for (int m_len = 0; m_len < three.Length; m_len++)
                            {
                                sw.Write(three_i[m_len] + ";");
                            }
                            sw.WriteLine();
                            for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++) 
                            {
                                sw.Write(bMax[bmax_len] + ";");
                                for (int m_len = 0; m_len < three.Length; m_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                 path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                 path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                 order ? aMax[amax_len] : bMax[bmax_len], 
                                                 order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void BaseMaxModGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo amax_di, one_dir;
            if (!order)
            {
                var x = aMax;
                aMax = bMax;
                bMax = x;
            }
            for (int f_len = 0; f_len < one.Length; f_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(one_i[f_len], one_i.Max()));
                for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                {
                    amax_di = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                    for (int m_len = 0; m_len < three.Length; m_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + amax_di.Name + "\\" + Choice[2] + "_" + AddZeros(three_i[m_len], three_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write((order ? "b_max" : "a_max") + "\\" + Two.Split('\\')[0] + ";");
 
                            for (int d_len = 0; d_len < two.Length; d_len++)
                            {
                                sw.Write(two_i[d_len] + ";");
                            }
                            sw.WriteLine();
                            for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++) 
                            {
                                sw.Write(bMax[bmax_len] + ";");
                                for (int d_len = 0; d_len < two.Length; d_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                  path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                  path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                  order ? aMax[amax_len] : bMax[bmax_len], 
                                                  order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }

        protected virtual void ExpBaseModGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            DirectoryInfo one_dir, two_dir;
            if (!order)
            {
                var x = one;
                one = three;
                three = x;
                
                var y = One;
                One = Three;
                Three = y;
            }
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            for (int d_len = 0; d_len < two.Length; d_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(two_i[d_len], two_i.Max()));
                for (int f_len = 0; f_len < one.Length; f_len++)
                {
                    two_dir = Directory.CreateDirectory(di + "\\" + one_dir.Name +  "\\" + Choice[1] + "_" + AddZeros(one_i[f_len], one_i.Max()));
                    for  (int m_len = 0; m_len < three.Length; m_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + two_dir.Name + "\\" + Choice[2] + "_" + AddZeros(three_i[m_len], three_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write("b_max" + "\\" + "a_max;");
 
                            for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                            {
                                sw.Write(aMax[amax_len] + ";");
                            }
                            sw.WriteLine();
                            for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                            {
                                sw.Write(bMax[bmax_len] + ";");
                                for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                                {
                                    sw.Write(Watch(path + "\\" + 
                                                   (order ? One + "\\" + one_i[f_len] : Three + "\\" + three_i[m_len]) + ".txt",
                                                 path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                 path + "\\" + 
                                                 (order ? Three + "\\" + three_i[m_len] : One + "\\" + one_i[f_len]) + ".txt", 
                                                 aMax[amax_len], bMax[bmax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void ExpBaseMax(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo one_dir, two_dir;
            if (!order)
            {
                var x = aMax;
                aMax = bMax;
                bMax = x;
            }
            for (int d_len = 0; d_len < two.Length; d_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(two_i[d_len], two_i.Max()));
                for (int f_len = 0; f_len < one.Length; f_len++)
                {
                    two_dir = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(one_i[f_len], one_i.Max()));
                    for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + two_dir.Name + "\\" + Choice[2] + "_" + AddZeros(aMax[amax_len], aMax.Max()) + ".csv", FileMode.Create, FileAccess.Write);
                        
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write((order ? "b_max" : "a_max") + "\\" + Three.Split('\\')[0] + ";");

                            for (int m_len = 0; m_len < three.Length; m_len++)
                            {
                                sw.Write(three_i[m_len] + ";");
                            }
                            sw.WriteLine();
                            for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                            {
                                sw.Write(bMax[bmax_len] + ";");
                                for (int m_len = 0; m_len < three.Length; m_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                 path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                 path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                 order ? aMax[amax_len] : bMax[bmax_len], 
                                                 order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void ExpModMaxGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo one_dir, two_dir;
            if (!order)
            {
                var x = aMax;
                aMax = bMax;
                bMax = x;
            }
            for (int d_len = 0; d_len < two.Length; d_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(two_i[d_len], two_i.Max()));
                for (int m_len = 0; m_len < three.Length; m_len++)
                {
                    two_dir = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(three_i[m_len], three_i.Max()));
                    for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + two_dir.Name + "\\" + Choice[2] + "_" + AddZeros(aMax[amax_len], aMax.Max()) + ".csv", FileMode.Create, FileAccess.Write);
                        
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write((order ? "b_max" : "a_max") + "\\" + One.Split('\\')[0] + ";");
 
                            for (int f_len = 0; f_len < one.Length; f_len++)
                            {
                                sw.Write(one_i[f_len] + ";");
                            }
                            sw.WriteLine();
                            for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                            {
                                sw.Write(bMax[bmax_len] + ";");
                                for (int f_len = 0; f_len < one.Length; f_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                  path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                  path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                  order ? aMax[amax_len] : bMax[bmax_len], 
                                                  order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void ExpMaxBaseGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo amax_di, one_dir;
            if (!order)
            {
               var x = aMax;
                aMax = bMax;
                bMax = x;
            }
            for (int d_len = 0; d_len < two.Length; d_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(two_i[d_len], two_i.Max()) );
                for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                {
                    amax_di = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                    for (int f_len = 0; f_len < one.Length; f_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + amax_di.Name + "\\" + Choice[2] + "_" + AddZeros(one_i[f_len], one_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write((order ? "b_max" : "a_max") + "\\" + Three.Split('\\')[0] + ";");
                            
                            for (int m_len = 0; m_len < three.Length; m_len++)
                            {
                                sw.Write(three_i[m_len] + ";");
                            }
                            sw.WriteLine();
                            for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++) 
                            {
                                sw.Write(bMax[bmax_len] + ";");
                                for (int m_len = 0; m_len < three.Length; m_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                 path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                 path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                 order ? aMax[amax_len] : bMax[bmax_len], 
                                                 order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void ExpMaxModGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo amax_di, one_dir;
            if (!order)
            {
                var x = aMax;
                aMax = bMax;
                bMax = x;
            }
            for (int d_len = 0; d_len < two.Length; d_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(two_i[d_len], two_i.Max()));
                for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                {
                    amax_di = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                    for (int m_len = 0; m_len < three.Length; m_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + amax_di.Name + "\\" + Choice[2] + "_" + AddZeros(three_i[m_len], three_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write((order ? "b_max" : "a_max") + "\\" + One.Split('\\')[0] + ";");
 
                            for (int f_len = 0; f_len < one.Length; f_len++)
                            {
                                sw.Write(one_i[f_len] + ";");
                            }
                            sw.WriteLine();
                            for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++) 
                            {
                                sw.Write(bMax[bmax_len] + ";");
                                for (int f_len = 0; f_len < one.Length; f_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                  path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                  path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                  order ? aMax[amax_len] : bMax[bmax_len], 
                                                  order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void ExpMaxMaxGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo amax_di, one_dir;
            if (!order)
            {
                var x = aMax;
                aMax = bMax;
                bMax = x; 
            }
            for (int d_len = 0; d_len < two.Length; d_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(two_i[d_len], two_i.Max()));
                for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                {
                    amax_di = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                    for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + amax_di.Name + "\\" + Choice[2] + "_" + AddZeros(bMax[bmax_len], bMax.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write(One.Split('\\')[0] + "\\" + Three.Split('\\')[0] + ";");
                            
                            for (int m_len = 0; m_len < three.Length; m_len++)
                            {
                                sw.Write(three_i[m_len] + ";");
                            }
                            sw.WriteLine();
                            for (int f_len = 0; f_len < one.Length; f_len++)
                            {
                                sw.Write(one_i[f_len] + ";");
                                for (int m_len = 0; m_len < three.Length; m_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                 path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                 path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                 order ? aMax[amax_len] : bMax[bmax_len], 
                                                 order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void ModBaseExpGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            DirectoryInfo one_dir, two_dir;
            if (!order)
            {
                var x = one;
                one = two;
                two = x;
                
                var y = One;
                One = Two;
                Two = y;
            }
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            for  (int m_len = 0; m_len < three.Length; m_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(three_i[m_len], three_i.Max()));
                for (int f_len = 0; f_len < one.Length; f_len++)
                {
                    two_dir = Directory.CreateDirectory(di + "\\" + one_dir.Name +  "\\" + Choice[1] + "_" + AddZeros(one_i[f_len], one_i.Max()));
                    for (int d_len = 0; d_len < two.Length; d_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + two_dir.Name + "\\" + Choice[2] + "_" + AddZeros(two_i[d_len], two_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write("b_max" + "\\" + "a_max;");
 
                            for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                            {
                                sw.Write(aMax[amax_len] + ";");
                            }
                            sw.WriteLine();
                            for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                            {
                                sw.Write(bMax[bmax_len] + ";");
                                for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                                {
                                    sw.Write(Watch(path + "\\" + 
                                                   (order ? One + "\\" + one_i[f_len] : Two + "\\" + two_i[d_len]) + ".txt",
                                                 path + "\\" + 
                                                 (order ? Two + "\\" + two_i[d_len] : One + "\\" + one_i[f_len]) + ".txt",
                                                 path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                 aMax[amax_len], bMax[bmax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void ModBaseMaxGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo one_dir, two_dir;
            if (!order)
            {
                var x = aMax;
                aMax = bMax;
                bMax = x;
            }
            for (int m_len = 0; m_len < three.Length; m_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(three_i[m_len], three_i.Max()));
                for (int f_len = 0; f_len < one.Length; f_len++)
                {
                    two_dir = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(one_i[f_len], one_i.Max()));
                    for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + two_dir.Name + "\\" + Choice[2] + "_" + AddZeros(aMax[amax_len], aMax.Max()) + ".csv", FileMode.Create, FileAccess.Write);
                        
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write((order ? "b_max" : "a_max") + "\\" + Two.Split('\\')[0] + ";");
 
                            for (int d_len = 0; d_len < two.Length; d_len++)
                            {
                                sw.Write(two_i[d_len] + ";");
                            }
                            sw.WriteLine();
                            for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                            {
                                sw.Write(bMax[bmax_len] + ";");
                                for (int d_len = 0; d_len < two.Length; d_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                  path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                  path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                  order ? aMax[amax_len] : bMax[bmax_len], 
                                                  order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void ModExpMaxGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo one_dir, two_dir;
            if (!order)
            {
                var x = aMax;
                aMax = bMax;
                bMax = x;
            }
            for (int m_len = 0; m_len < three.Length; m_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(three_i[m_len], three_i.Max()));
                for (int d_len = 0; d_len < two.Length; d_len++)
                {
                    two_dir = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(two_i[d_len], two_i.Max()));
                    for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + two_dir.Name + "\\" + Choice[2] + "_" + AddZeros(aMax[amax_len], aMax.Max()) + ".csv", FileMode.Create, FileAccess.Write);
                        
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write((order ? "b_max" : "a_max") + "\\" + One.Split('\\')[0] + ";");
 
                            for (int f_len = 0; f_len < one.Length; f_len++)
                            {
                                sw.Write(one_i[f_len] + ";");
                            }
                            sw.WriteLine();
                            for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                            {
                                sw.Write(bMax[bmax_len] + ";");
                                for (int f_len = 0; f_len < one.Length; f_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                  path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                  path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                  order ? aMax[amax_len] : bMax[bmax_len], 
                                                  order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void ModMaxBaseGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo amax_di, one_dir;
            if (!order)
            {
               var x = aMax;
                aMax = bMax;
                bMax = x;
            }
            for (int m_len = 0; m_len < three.Length; m_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(three_i[m_len], three_i.Max()) );
                for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                {
                    amax_di = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                    for (int f_len = 0; f_len < one.Length; f_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + amax_di.Name + "\\" + Choice[2] + "_" + AddZeros(one_i[f_len], one_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write((order ? "b_max" : "a_max") + "\\" + Two.Split('\\')[0] + ";");
                            
                            for (int d_len = 0; d_len < two.Length; d_len++)
                            {
                                sw.Write(two_i[d_len] + ";");
                            }
                            sw.WriteLine();
                            for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++) 
                            {
                                sw.Write(bMax[bmax_len] + ";");
                                for (int d_len = 0; d_len < two.Length; d_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                 path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                 path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                 order ? aMax[amax_len] : bMax[bmax_len], 
                                                 order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void ModMaxExpGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo amax_di, one_dir;
            if (!order)
            {
                var x = aMax;
                aMax = bMax;
                bMax = x;
            }
            for (int m_len = 0; m_len < three.Length; m_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(three_i[m_len], three_i.Max()) );
                for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                {
                    amax_di = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                    for (int d_len = 0; d_len < two.Length; d_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + amax_di.Name + "\\" + Choice[2] + "_" + AddZeros(two_i[d_len], two_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write((order ? "b_max" : "a_max") + "\\" + One.Split('\\')[0] + ";");
 
                            for (int f_len = 0; f_len < one.Length; f_len++)
                            {
                                sw.Write(two_i[d_len] + ";");
                            }
                            sw.WriteLine();
                            for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++) 
                            {
                                sw.Write(bMax[bmax_len] + ";");
                                for (int f_len = 0; f_len < one.Length; f_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                  path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                  path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                  order ? aMax[amax_len] : bMax[bmax_len], 
                                                  order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void ModMaxMaxGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo amax_di, one_dir;
            if (!order)
            {
                var x = aMax;
                aMax = bMax;
                bMax = x; 
            }
            for (int m_len = 0; m_len < three.Length; m_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(three_i[m_len], three_i.Max()));
                for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
                {
                    amax_di = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                    for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + amax_di.Name + "\\" + Choice[2] + "_" + AddZeros(bMax[bmax_len], bMax.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write(One.Split('\\')[0] + "\\" + Two.Split('\\')[0] + ";");
                            
                            for (int d_len = 0; d_len < two.Length; d_len++)
                            {
                                sw.Write(two_i[d_len] + ";");
                            }
                            sw.WriteLine();
                            for (int f_len = 0; f_len < one.Length; f_len++)
                            {
                                sw.Write(one_i[f_len] + ";");
                                for (int d_len = 0; d_len < two.Length; d_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                 path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                 path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                 order ? aMax[amax_len] : bMax[bmax_len], 
                                                 order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void MaxBaseModGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo amax_di, bmax_di;
            if (!order)
            {
                var x = aMax;
                aMax = bMax;
                bMax = x;
            }
            for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
            {
                amax_di = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                for (int f_len = 0; f_len < one.Length; f_len++)
                {
                    bmax_di = Directory.CreateDirectory(di + "\\" + amax_di.Name + "\\" + Choice[1] + "_" + AddZeros(one_i[f_len], one_i.Max()));
                    for (int m_len = 0; m_len < three.Length; m_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + amax_di.Name + "\\" + bmax_di.Name + "\\" + Choice[2] + "_" + AddZeros(three_i[m_len], three_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
                        
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write(Two.Split('\\')[0] + "\\" + (order ? "b_max" : "a_max") + ";");
                            
                            for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++) 
                            {
                                sw.Write(bMax[bmax_len] + ";");
                            }
                            sw.WriteLine();
                            for (int d_len = 0; d_len < two.Length; d_len++)
                            {
                                sw.Write(two_i[d_len] + ";");
                                for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++) 
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                 path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                 path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                 order ? aMax[amax_len] : bMax[bmax_len], 
                                                 order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void MaxExpBaseGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo amax_di, one_dir;
            if (!order)
            {
               var x = aMax;
                aMax = bMax;
                bMax = x;
            }
            for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                for (int d_len = 0; d_len < two.Length; d_len++)
                {
                    amax_di = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(two_i[d_len], two_i.Max()));
                    for (int f_len = 0; f_len < one.Length; f_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + amax_di.Name + "\\" + Choice[2] + "_" + AddZeros(one_i[f_len], one_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write((order ? "b_max" : "a_max") + "\\" + Three.Split('\\')[0] + ";");
                            
                            for (int m_len = 0; m_len < three.Length; m_len++)
                            {
                                sw.Write(three_i[m_len] + ";");
                            }
                            sw.WriteLine();
                            for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++) 
                            {
                                sw.Write(bMax[bmax_len] + ";");
                                for (int m_len = 0; m_len < three.Length; m_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                 path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                 path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                 order ? aMax[amax_len] : bMax[bmax_len], 
                                                 order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void MaxExpModGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo amax_di, one_dir;
            if (!order)
            {
                var x = aMax;
                aMax = bMax;
                bMax = x;
            }
            for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                for (int d_len = 0; d_len < two.Length; d_len++)
                {
                    amax_di = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(two_i[d_len], two_i.Max()));
                    for (int m_len = 0; m_len < three.Length; m_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + amax_di.Name + "\\" + Choice[2] + "_" + AddZeros(three_i[m_len], three_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write((order ? "b_max" : "a_max") + "\\" + One.Split('\\')[0] + ";");
 
                            for (int f_len = 0; f_len < one.Length; f_len++)
                            {
                                sw.Write(one_i[f_len] + ";");
                            }
                            sw.WriteLine();
                            for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++) 
                            {
                                sw.Write(bMax[bmax_len] + ";");
                                for (int f_len = 0; f_len < one.Length; f_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                  path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                  path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                  order ? aMax[amax_len] : bMax[bmax_len], 
                                                  order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void MaxModBaseGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo amax_di, bmax_di;
            if (!order)
            {
                var x = aMax;
                aMax = bMax;
                bMax = x;
            }
            for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
            {
                amax_di = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                for (int m_len = 0; m_len < three.Length; m_len++)
                {
                    bmax_di = Directory.CreateDirectory(di + "\\" + amax_di.Name + "\\" + Choice[1] + "_" + AddZeros(three_i[m_len], three_i.Max()));
                    for (int f_len = 0; f_len < one.Length; f_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + amax_di.Name + "\\" + bmax_di.Name + "\\" + Choice[2] + "_" + AddZeros(one_i[f_len], one_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
                        
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write(Two.Split('\\')[0] + "\\" + (order ? "b_max" : "a_max") + ";");
                            
                            for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++) 
                            {
                                sw.Write(bMax[bmax_len] + ";");
                            }
                            sw.WriteLine();
                            for (int d_len = 0; d_len < two.Length; d_len++)
                            {
                                sw.Write(two_i[d_len] + ";");
                                for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++) 
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                 path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                 path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                 order ? aMax[amax_len] : bMax[bmax_len], 
                                                 order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void MaxModExpGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo amax_di, one_dir;
            if (!order)
            {
                var x = aMax;
                aMax = bMax;
                bMax = x;
            }
            for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                for (int m_len = 0; m_len < three.Length; m_len++)
                {
                    amax_di = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(three_i[m_len], three_i.Max()));
                    for (int d_len = 0; d_len < two.Length; d_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + amax_di.Name + "\\" + Choice[2] + "_" + AddZeros(two_i[d_len], two_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write((order ? "b_max" : "a_max") + "\\" + One.Split('\\')[0] + ";");
 
                            for (int f_len = 0; f_len < one.Length; f_len++)
                            {
                                sw.Write(one_i[f_len] + ";");
                            }
                            sw.WriteLine();
                            for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++) 
                            {
                                sw.Write(bMax[bmax_len] + ";");
                                for (int f_len = 0; f_len < one.Length; f_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                  path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                  path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                  order ? aMax[amax_len] : bMax[bmax_len], 
                                                  order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void MaxModMaxGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo amax_di, one_dir;
            if (!order)
            {
                var x = aMax;
                aMax = bMax;
                bMax = x; 
            }
            for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                for (int m_len = 0; m_len < three.Length; m_len++)
                {
                    amax_di = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(three_i[m_len], three_i.Max()));
                    for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + amax_di.Name + "\\" + Choice[2] + "_" + AddZeros(bMax[bmax_len], bMax.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write(One.Split('\\')[0] + "\\" + Two.Split('\\')[0] + ";");
                            
                            for (int d_len = 0; d_len < two.Length; d_len++)
                            {
                                sw.Write(two_i[d_len] + ";");
                            }
                            sw.WriteLine();
                            for (int f_len = 0; f_len < one.Length; f_len++)
                            {
                                sw.Write(one_i[f_len] + ";");
                                for (int d_len = 0; d_len < two.Length; d_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                 path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                 path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                 order ? aMax[amax_len] : bMax[bmax_len], 
                                                 order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void MaxMaxBaseGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo amax_di, one_dir;
            if (!order)
            {
                var x = aMax;
                aMax = bMax;
                bMax = x; 
            }
            for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                {
                    amax_di = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(bMax[bmax_len], bMax.Max()));
                    for (int f_len = 0; f_len < one.Length; f_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + amax_di.Name + "\\" + Choice[2] + "_" + AddZeros(one_i[f_len], one_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write(Two.Split('\\')[0] + "\\" + Three.Split('\\')[0] + ";");
                            
                            for (int m_len = 0; m_len < three.Length; m_len++)
                            {
                                sw.Write(three_i[m_len] + ";");
                            }
                            sw.WriteLine();
                            for (int d_len = 0; d_len < two.Length; d_len++)
                            {
                                sw.Write(two_i[d_len] + ";");
                                for (int m_len = 0; m_len < three.Length; m_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                 path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                 path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                 order ? aMax[amax_len] : bMax[bmax_len], 
                                                 order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        protected virtual void MaxMaxExpGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
         {
             string path = Directory.GetCurrentDirectory();
             int[] one_i = MakeDigits(one);
             int[] two_i = MakeDigits(two);
             int[] three_i = MakeDigits(three);
             DirectoryInfo amax_di, bmax_di;
             if (!order)
             {
                 var x = aMax;
                 aMax = bMax;
                 bMax = x;
             }
             for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
             {
                 amax_di = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                 for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                 {
                     bmax_di = Directory.CreateDirectory(di + "\\" + amax_di.Name + "\\" + Choice[1] + "_" + AddZeros(bMax[bmax_len], bMax.Max()));
                     for (int d_len = 0; d_len < two.Length; d_len++)
                     {
                         FileStream fin = new FileStream(di + "\\" + amax_di.Name + "\\" + bmax_di.Name + "\\" + Choice[2] + "_" + AddZeros(two_i[d_len], two_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                         using (StreamWriter sw = new StreamWriter(fin))
                         {
                             sw.Write(One.Split('\\')[0] + "\\" + Three.Split('\\')[0] + ";");
 
                             for (int m_len = 0; m_len < three.Length; m_len++)
                             {
                                 sw.Write(three_i[m_len] + ";");
                             }
                             sw.WriteLine();
                             for (int f_len = 0; f_len < one.Length; f_len++)
                             {
                                 sw.Write(one_i[f_len] + ";");
                                 for (int m_len = 0; m_len < three.Length; m_len++)
                                 {
                                     sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                  path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                  path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                  order ? aMax[amax_len] : bMax[bmax_len], 
                                                  order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                 }
                                 sw.WriteLine();
                             }
                         }
                         fin.Close();
                     }
                 }
             }
         }
        
        protected virtual void MaxMaxModGen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, bool order = true)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            DirectoryInfo amax_di, one_dir;
            if (!order)
            {
                var x = aMax;
                aMax = bMax;
                bMax = x; 
            }
            for (int amax_len = 0; amax_len < aMax.Count; amax_len++)
            {
                one_dir = Directory.CreateDirectory(di + "\\" + Choice[0] + "_" + AddZeros(aMax[amax_len], aMax.Max()));
                for (int bmax_len = 0; bmax_len < bMax.Count; bmax_len++)
                {
                    amax_di = Directory.CreateDirectory(di + "\\" + one_dir.Name + "\\" + Choice[1] + "_" + AddZeros(bMax[bmax_len], bMax.Max()));
                    for (int m_len = 0; m_len < three.Length; m_len++)
                    {
                        FileStream fin = new FileStream(di + "\\" + one_dir.Name + "\\" + amax_di.Name + "\\" + Choice[2] + "_" + AddZeros(three_i[m_len], three_i.Max()) + ".csv", FileMode.Create, FileAccess.Write);
 
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write(One.Split('\\')[0] + "\\" + Two.Split('\\')[0] + ";");
                            
                            for (int d_len = 0; d_len < two.Length; d_len++)
                            {
                                sw.Write(two_i[d_len] + ";");
                            }
                            sw.WriteLine();
                            for (int f_len = 0; f_len < one.Length; f_len++)
                            {
                                sw.Write(one_i[f_len] + ";");
                                for (int d_len = 0; d_len < two.Length; d_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                 path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                 path + "\\" + Three + "\\" + three_i[m_len] + ".txt", 
                                                 order ? aMax[amax_len] : bMax[bmax_len], 
                                                 order ? bMax[bmax_len] : aMax[amax_len]) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
            }
        }
        
        #endregion
 
        public virtual void Create_Result()
        {
            SendStatus();
            string path = Directory.GetCurrentDirectory();
            string[] founds = Directory.GetFiles(path + "\\Base\\" + Found);
            string[] degrees = Directory.GetFiles(path + "\\Exponent\\" + Degree);
            string[] mods = Directory.GetFiles(path + "\\Modulus\\" + Mod);
 
            DirectoryInfo di = CD(Choice, Found, Degree, Mod);
            var di_short = "Results\\" + di.Name;
            switch (Choice[0])
            {
                case "Base":
                    switch (Choice[1])
                    {
                        case "Exponent":
                            switch (Choice[2])
                            {
                                case "Modulus":
                                    BaseExpModGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "aMax":
                                    BaseExpMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "bMax":
                                    BaseExpMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                            }
                            break;
                        case "Modulus":
                            switch (Choice[2])
                            {
                                case "Exponent":
                                    BaseExpModGen(founds, mods, degrees , di_short, "Base\\" + Found, "Modulus\\" + Mod, "Exponent\\" + Degree, false);
                                    break;
                                case "aMax":
                                    BaseModMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "bMax":
                                    BaseModMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                            }
                            break;
                        case "aMax":
                            switch (Choice[2])
                            {
                                case "Exponent":
                                    BaseMaxExpGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "Modulus":
                                    BaseMaxModGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "bMax":
                                    BaseMaxMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                            }
                            break;
                        case "bMax":
                            switch (Choice[2])
                            {
                                case "Exponent":
                                    BaseMaxExpGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                                case "Modulus":
                                    BaseMaxModGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                                case "aMax":
                                    BaseMaxMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                            }
                            break;
                    }
                    break;
                case "Exponent":
                    switch (Choice[1])
                    {
                        case "Base":
                            switch (Choice[2])
                            {
                                case "Modulus":
                                    ExpBaseModGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "aMax":
                                    ExpBaseMax(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "bMax":
                                    ExpBaseMax(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                            }
                            break;
                        case "Modulus":
                            switch (Choice[2])
                            {
                                case "Base":
                                    ExpBaseModGen(mods, degrees, founds, di_short, "Modulus\\" + Mod, "Exponent\\" + Degree, "Base\\" + Found, false);
                                    break;
                                case "aMax":
                                    ExpModMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "bMax":
                                    ExpModMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                            }
                            break;
                        case "aMax":
                            switch (Choice[2])
                            {
                                case "Base":
                                    ExpMaxBaseGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "Modulus":
                                    ExpMaxModGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "bMax":
                                    ExpMaxMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                            }
                            break;
                        case "bMax":
                            switch (Choice[2])
                            {
                                case "Base":
                                    ExpMaxBaseGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                                case "Modulus":
                                    ExpMaxModGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                                case "aMax":
                                    ExpMaxMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                            }
                            break;
                    }
                    break;
                case "Modulus":
                    switch (Choice[1])
                    {
                        case "Base":
                            switch (Choice[2])
                            {
                                case "Exponent":
                                    ModBaseExpGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "aMax":
                                    ModBaseMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "bMax":
                                    ModBaseMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                            }
                            break;
                        case "Exponent":
                            switch (Choice[2])
                            {
                                case "Base":
                                    ModBaseExpGen(degrees, founds, mods, di_short, "Exponent\\" + Degree, "Base\\" + Found, "Modulus\\" + Mod, false);
                                    break;
                                case "aMax":
                                    ModExpMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "bMax":
                                    ModExpMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                            }
                            break;
                        case "aMax":
                            switch (Choice[2])
                            {
                                case "Base":
                                    ModMaxBaseGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "Exponent":
                                    ModMaxExpGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "bMax":
                                    ModMaxMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                            }
                            break;
                        case "bMax":
                            switch (Choice[2])
                            {
                                case "Base":
                                    ModMaxBaseGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                                case "Exponent":
                                    ModMaxExpGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                                case "aMax":
                                    ModMaxMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                            }
                            break;
                    }
                    break;
                case "aMax":
                    switch (Choice[1])
                    {
                        case "Base":
                            switch (Choice[2])
                            {
                                case "Exponent":
                                    MaxBaseExpGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "Modulus":
                                    MaxBaseModGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "bMax":
                                    MaxBaseMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                            }
                            break;
                        case "Exponent":
                            switch (Choice[2])
                            {
                                case "Base":
                                    MaxExpBaseGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "Modulus":
                                    MaxExpModGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "bMax":
                                    MaxExpMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                            }
                            break;
                        case "Modulus":
                            switch (Choice[2])
                            {
                                case "Base":
                                    MaxModBaseGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "Exponent":
                                    MaxModExpGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "bMax":
                                    MaxModMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                            }
                            break;
                        case "bMax":
                            switch (Choice[2])
                            {
                                case "Base":
                                    MaxMaxBaseGen(founds, degrees, mods, di_short, "Base\\" + Found,  "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "Exponent":
                                    MaxMaxExpGen(founds, degrees, mods, di_short, "Base\\" + Found,  "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                                case "Modulus":
                                    MaxMaxModGen(founds, degrees, mods, di_short, "Base\\" + Found,  "Exponent\\" + Degree, "Modulus\\" + Mod);
                                    break;
                            }
                            break;
                    }
                    break;
                case "bMax":
                    switch (Choice[1])
                    {
                        case "Base":
                            switch (Choice[2])
                            {
                                case "Exponent":
                                    MaxBaseExpGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                                case "Modulus":
                                    MaxBaseModGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                                case "aMax":
                                    MaxBaseMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                            }
                            break;
                        case "Exponent":
                            switch (Choice[2])
                            {
                                case "Base":
                                    MaxExpBaseGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                                case "Modulus":
                                    MaxExpModGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                                case "aMax":
                                    MaxExpMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                            }
                            break;
                        case "Modulus":
                            switch (Choice[2])
                            {
                                case "Base":
                                    MaxModBaseGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                                case "Exponent":
                                    MaxModExpGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                                case "aMax":
                                    MaxModMaxGen(founds, degrees, mods, di_short, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                            }
                            break;
                        case "aMax":
                            switch (Choice[2])
                            {
                                case "Base":
                                    MaxMaxBaseGen(founds, degrees, mods, di_short, "Base\\" + Found,  "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                                case "Exponent":
                                    MaxMaxExpGen(founds, degrees, mods, di_short, "Base\\" + Found,  "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                                case "Modulus":
                                    MaxMaxModGen(founds, degrees, mods, di_short, "Base\\" + Found,  "Exponent\\" + Degree, "Modulus\\" + Mod, false);
                                    break;
                            }
                            break;
                    }
                    break;
            }
            SendStatus();
        }

        protected double Watch(string Found_file, string Degree_file, string Mod_file, BigInteger amax, BigInteger bmax)
        {
            Stopwatch stw = new Stopwatch();
            List<double> sums = new List<double>();
            double f_sum = 0;
            int f_count = 0;

            using (StreamReader srf = new StreamReader(Found_file))
            {
                string found;
                while ((found = srf.ReadLine()) != null)
                {
                    double d_sum = 0;
                    int d_count = 0;
                    using (StreamReader srd = new StreamReader(Degree_file))
                    {
                        string degree;
                        while ((degree = srd.ReadLine()) != null)
                        {
                            double m_sum = 0;
                            int m_count = 0;
                            using (StreamReader srm = new StreamReader(Mod_file))
                            {
                                string mod;
                                while ((mod = srm.ReadLine()) != null)
                                {
                                    stw.Start();
                                    LoopFunc(BigInteger.Parse(found), BigInteger.Parse(degree), BigInteger.Parse(mod), amax, bmax);
                                    stw.Stop();

                                    m_sum += stw.Elapsed.TotalMilliseconds;
                                    m_count++;
                                    stw.Reset();
                                }

                            }
                            d_sum += m_sum / m_count;
                            d_count++;
                        }
                    }
                    f_sum += d_sum / d_count;
                    f_count++;
                }
            }
            return f_sum / f_count;
        }
    }

    abstract class WindowPow 
    {
        public Inverse InvMethod { get; set; }
        public Multiply MultMethod { get; set; }

        public string Found { get; set; }
        public string Degree { get; set; }
        public string Mod { get; set; }
        public string Choice { get; set; }
        public string Window { get; set; }
        public bool Table { get; set; }

        public WindowPow(string found, string degree, string mod, string by, string window, bool table)
        {
            Found = found;
            Degree = degree;
            Mod = mod;
            Choice = by;
            Window = window;
            Table = table;
        }

        public abstract string Name();

        protected abstract void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime);

        protected void SendStatus()
        {
            (Application.OpenForms[0] as MainForm).ChangeStatus(Name());
        }

        private int[] MakeDigits(string[] mas)
        {
            int[] mas_i = new int[mas.Length]; ;

            for (int i = 0; i < mas.Length; i++)
            {
                var s = mas[i].Split('\\');
                mas_i[i] = Int32.Parse(s[s.Length - 1].Split('.')[0]);
            }
            Array.Sort(mas_i);
            return mas_i;
        }

        private void Gen(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, string w, bool Table)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            var w_i = GenFunctions.ReadString(w);
            DirectoryInfo o_di;
            for (int f_len = 0; f_len < one.Length; f_len++)
            {
                o_di = Directory.CreateDirectory(di + "\\" + One.Split('\\')[0] + "_" + one_i[f_len]);
                double table_time;
                /*
                if (Table)
                {
                    for (int d_len = 0; d_len < two.Length; d_len++)
                    {
                        FileStream fin = new FileStream(o_di.FullName + "\\" + Two.Split('\\')[0] + "_" + two_i[d_len] + ".csv", FileMode.Create, FileAccess.Write);

                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write("Win" + "\\" + Three.Split('\\')[0] + ";");

                            for (int m_len = 0; m_len < three.Length; m_len++)
                            {
                                sw.Write(three_i[m_len] + ";");
                            }
                            sw.WriteLine();
                            for (int w_len = 0; w_len < w_i.Count; w_len++)
                            {
                                sw.Write(w_i[w_len] + ";");
                                for (int m_len = 0; m_len < three.Length; m_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                   path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                   path + "\\" + Three + "\\" + three_i[m_len] + ".txt", w_i[w_len], out table_time) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
                else
                {*/
                string str = "";

                if (Table)
                    str = "_full_time";

                for (int d_len = 0; d_len < two.Length; d_len++)
                    {
                            FileStream fin = new FileStream(o_di.FullName + "\\" + Two.Split('\\')[0] + "_" + two_i[d_len] + str + ".csv", FileMode.Create, FileAccess.Write);
                    FileStream table = new FileStream(o_di.FullName + "\\" + Two.Split('\\')[0] + "_" + two_i[d_len] + "_table_time" + ".csv", FileMode.Create, FileAccess.Write);
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            using (StreamWriter sw2 = new StreamWriter(table))
                            {
                                sw.Write("Win" + "\\" + Three.Split('\\')[0] + ";");
                                sw2.Write("Win" + "\\" + Three.Split('\\')[0] + ";");
                                for (int m_len = 0; m_len < three.Length; m_len++)
                                {
                                    sw.Write(three_i[m_len] + ";");
                                    sw2.Write(three_i[m_len] + ";");
                                }
                                sw.WriteLine();
                                sw2.WriteLine();
                                for (int w_len = 0; w_len < w_i.Count; w_len++)
                                {
                                    sw.Write(w_i[w_len] + ";");
                                    sw2.Write(w_i[w_len] + ";");
                                    for (int m_len = 0; m_len < three.Length; m_len++)
                                    {
                                        if (Table)
                                            sw.Write((Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                       path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                       path + "\\" + Three + "\\" + three_i[m_len] + ".txt", w_i[w_len], out table_time)) + ";");
                                        else
                                            sw.Write((Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                   path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                   path + "\\" + Three + "\\" + three_i[m_len] + ".txt", w_i[w_len], out table_time) - table_time) + ";");
                                    sw2.Write(table_time + ";");
                                    }
                                    sw.WriteLine();
                                    sw2.WriteLine();
                                }
                            }
                        }
                        fin.Close();
                        table.Close();
                    }
                //}
            }
        }

        private void GenWin(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, string w, bool Table)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            var w_i = GenFunctions.ReadString(w);
            DirectoryInfo w_di;
            for (int w_len = 0; w_len < w_i.Count; w_len++)
            {
                w_di = Directory.CreateDirectory(di + "\\Window" + w_i[w_len]);
                double table_time;
                /*
                if (Table)
                {
                    for (int f_len = 0; f_len < one.Length; f_len++)
                    {
                        FileStream fin = new FileStream(w_di.FullName + "\\" + One.Split('\\')[0] + "_" + one_i[f_len] + ".csv", FileMode.Create, FileAccess.Write);

                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write(Two.Split('\\')[0] + "\\" + Three.Split('\\')[0] + ";");

                            for (int m_len = 0; m_len < three.Length; m_len++)
                            {
                                sw.Write(three_i[m_len] + ";");
                            }
                            sw.WriteLine();
                            for (int d_len = 0; d_len < two.Length; d_len++)
                            {
                                sw.Write(two_i[d_len] + ";");
                                for (int m_len = 0; m_len < three.Length; m_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                   path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                   path + "\\" + Three + "\\" + three_i[m_len] + ".txt", w_i[w_len], out table_time) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
                else
                {*/
                string str = "";

                if (Table)
                    str = "_full_time";

                for (int f_len = 0; f_len < one.Length; f_len++)
                    {
                        FileStream fin = new FileStream(w_di.FullName + "\\" + One.Split('\\')[0] + "_" + one_i[f_len] + str + ".csv", FileMode.Create, FileAccess.Write);
                        FileStream table = new FileStream(w_di.FullName + "\\" + One.Split('\\')[0] + "_" + one_i[f_len] + "_table_time" + ".csv", FileMode.Create, FileAccess.Write);
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            using (StreamWriter sw2 = new StreamWriter(table))
                            {
                                sw.Write(Two.Split('\\')[0] + "\\" + Three.Split('\\')[0] + ";");
                                sw2.Write(Two.Split('\\')[0] + "\\" + Three.Split('\\')[0] + ";");
                                for (int m_len = 0; m_len < three.Length; m_len++)
                                {
                                    sw.Write(three_i[m_len] + ";");
                                    sw2.Write(three_i[m_len] + ";");
                                }
                                sw.WriteLine();
                                sw2.WriteLine();
                                for (int d_len = 0; d_len < two.Length; d_len++)
                                {
                                    sw.Write(two_i[d_len] + ";");
                                    sw2.Write(two_i[d_len] + ";");
                                    for (int m_len = 0; m_len < three.Length; m_len++)
                                    {
                                        if (Table)
                                        sw.Write((Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                       path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                       path + "\\" + Three + "\\" + three_i[m_len] + ".txt", w_i[w_len], out table_time)) + ";");
                                        else
                                        sw.Write((Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                       path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                       path + "\\" + Three + "\\" + three_i[m_len] + ".txt", w_i[w_len], out table_time) - table_time) + ";");

                                    sw2.Write(table_time + ";");
                                    }
                                    sw.WriteLine();
                                    sw2.WriteLine();
                                }
                            }
                        }
                        fin.Close();
                        table.Close();
                    }
                //}
            }
        }

        private void GenWin2(string[] one, string[] two, string[] three, string di, string One, string Two, string Three, string w, bool Table)
        {
            string path = Directory.GetCurrentDirectory();
            int[] one_i = MakeDigits(one);
            int[] two_i = MakeDigits(two);
            int[] three_i = MakeDigits(three);
            var w_i = GenFunctions.ReadString(w);
            DirectoryInfo o_di;
            for (int f_len = 0; f_len < one.Length; f_len++)
            {
                o_di = Directory.CreateDirectory(di + "\\" + One.Split('\\')[0] + "_" + one_i[f_len]);
                double table_time;
                /*
                if (Table)
                {
                    for (int w_len = 0; w_len < w_i.Count; w_len++)
                    {
                        FileStream fin = new FileStream(o_di.FullName + "\\" + One.Split('\\')[0] + "Window_" + w_i[w_len] + ".csv", FileMode.Create, FileAccess.Write);

                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            sw.Write(Two.Split('\\')[0] + "\\" + Three.Split('\\')[0] + ";");

                            for (int m_len = 0; m_len < three.Length; m_len++)
                            {
                                sw.Write(three_i[m_len] + ";");
                            }
                            sw.WriteLine();
                            for (int d_len = 0; d_len < two.Length; d_len++)
                            {
                                sw.Write(two_i[d_len] + ";");
                                for (int m_len = 0; m_len < three.Length; m_len++)
                                {
                                    sw.Write(Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                   path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                   path + "\\" + Three + "\\" + three_i[m_len] + ".txt", w_i[w_len], out table_time) + ";");
                                }
                                sw.WriteLine();
                            }
                        }
                        fin.Close();
                    }
                }
                else
                {*/
                string str = "";

                if (Table)
                    str = "_full_time";

                for (int w_len = 0; w_len < w_i.Count; w_len++)
                    {
                        FileStream fin = new FileStream(o_di.FullName + "\\" + One.Split('\\')[0] + "Window_" + w_i[w_len] + str + ".csv", FileMode.Create, FileAccess.Write);
                        FileStream table = new FileStream(o_di.FullName + "\\" + One.Split('\\')[0] + "Window_" + w_i[w_len] + "_table_time" + ".csv", FileMode.Create, FileAccess.Write);
                        using (StreamWriter sw = new StreamWriter(fin))
                        {
                            using (StreamWriter sw2 = new StreamWriter(table))
                            {
                                sw.Write(Two.Split('\\')[0] + "\\" + Three.Split('\\')[0] + ";");
                                sw2.Write(Two.Split('\\')[0] + "\\" + Three.Split('\\')[0] + ";");
                                for (int m_len = 0; m_len < three.Length; m_len++)
                                {
                                    sw.Write(three_i[m_len] + ";");
                                    sw2.Write(three_i[m_len] + ";");
                                }
                                sw.WriteLine();
                                sw2.WriteLine();
                                for (int d_len = 0; d_len < two.Length; d_len++)
                                {
                                    sw.Write(two_i[d_len] + ";");
                                    sw2.Write(two_i[d_len] + ";");
                                    for (int m_len = 0; m_len < three.Length; m_len++)
                                    {
                                        if (Table)
                                            sw.Write((Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                       path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                       path + "\\" + Three + "\\" + three_i[m_len] + ".txt", w_i[w_len], out table_time)) + ";");
                                        else
                                            sw.Write((Watch(path + "\\" + One + "\\" + one_i[f_len] + ".txt",
                                                       path + "\\" + Two + "\\" + two_i[d_len] + ".txt",
                                                       path + "\\" + Three + "\\" + three_i[m_len] + ".txt", w_i[w_len], out table_time) - table_time) + ";");
                                        sw2.Write(table_time + ";");
                                    }
                                    sw.WriteLine();
                                    sw2.WriteLine();
                                }
                            }
                        }
                        fin.Close();
                        table.Close();
                    }
                //}
            }
        }

        private DirectoryInfo CD(string by, string Found, string Degree, string Mod, string window)
        {
            DirectoryInfo di;
            string path = Directory.GetCurrentDirectory();
            var lengts = GenFunctions.ReadString(window);
            di = Directory.CreateDirectory(String.Format("{0}\\{1}\\{9}{2}_{3}_{4}_{5}_Window_{6}-{7}(" + by.Split(' ')[0][0] + "_" + by.Split(' ')[1][0] + ")#{8}",
                                         path, "Results", Name(), Mod.Split('#')[0], Found.Split('#')[0], Degree.Split('#')[0], lengts[0], lengts[lengts.Count - 1],
                                         DateTime.Now.ToLocalTime().ToString().Replace(':', '-'), AdditionalParameters.montFlag ? "Mont_" : ""));
            return di;
        }

        public void Create_Result()
        {
            SendStatus();
            string path = Directory.GetCurrentDirectory();

            string[] founds = Directory.GetFiles(path + "\\Base\\" + Found);
            string[] degrees = Directory.GetFiles(path + "\\Exponent\\" + Degree);
            string[] mods = Directory.GetFiles(path + "\\Modulus\\" + Mod);
            DirectoryInfo di = CD(Choice, Found, Degree, Mod, Window);
            switch (Choice.Split(' ')[0])
            {
                case "Base":
                    switch (Choice.Split(' ')[1])
                    {
                        case "Window":
                            GenWin2(founds, degrees, mods, di.FullName, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, Window, Table);
                            break;
                        case "Exponent":
                            Gen(founds, degrees, mods, di.FullName, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, Window, Table);
                            break;
                        case "Modulus":
                            Gen(founds, mods, degrees, di.FullName, "Base\\" + Found, "Modulus\\" + Mod, "Exponent\\" + Degree, Window, Table);
                            break;
                        default: goto smth_go_wrong;
                    }
                    break;

                case "Exponent":
                    switch (Choice.Split(' ')[1])
                    {
                        case "Window":
                            GenWin2(degrees, founds, mods, di.FullName, "Exponent\\" + Degree, "Base\\" + Found, "Modulus\\" + Mod, Window, Table);
                            break;
                        case "Base":
                            Gen(degrees, founds, mods, di.FullName, "Exponent\\" + Degree, "Base\\" + Found, "Modulus\\" + Mod, Window, Table);
                            break;
                        case "Modulus":
                            Gen(degrees, mods, founds, di.FullName, "Exponent\\" + Degree, "Modulus\\" + Mod, "Base\\" + Found, Window, Table);
                            break;
                        default: goto smth_go_wrong;
                    }
                    break;

                case "Modulus":
                    switch (Choice.Split(' ')[1])
                    {
                        case "Window":
                            GenWin2(mods, degrees, founds, di.FullName, "Modulus\\" + Mod, "Exponent\\" + Degree, "Base\\" + Found, Window, Table);
                            break;
                        case "Base":
                            Gen(mods, founds, degrees, di.FullName, "Modulus\\" + Mod, "Base\\" + Found, "Exponent\\" + Degree, Window, Table);
                            break;
                        case "Exponent":
                            Gen(mods, degrees, founds, di.FullName, "Modulus\\" + Mod, "Exponent\\" + Degree, "Base\\" + Found, Window, Table);
                            break;
                        default: goto smth_go_wrong;
                    }

                    break;


                case "Window":
                    switch (Choice.Split(' ')[1])
                    {
                        case "Exponent":
                            GenWin(degrees, founds, mods, di.FullName, "Exponent\\" + Degree, "Base\\" + Found, "Modulus\\" + Mod, Window, Table);
                            break;
                        case "Base":
                            GenWin(founds, degrees, mods, di.FullName, "Base\\" + Found, "Exponent\\" + Degree, "Modulus\\" + Mod, Window, Table);
                            break;
                        case "Modulus":
                            GenWin(mods, degrees, founds, di.FullName, "Modulus\\" + Mod, "Exponent\\" + Degree, "Base\\" + Found, Window, Table);
                            break;
                        default: goto smth_go_wrong;
                    }
                    break;

                default:
                smth_go_wrong:
                    MessageBox.Show("Wrong parametrs", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                    break;
            }
            SendStatus();
        }

        public void Create_Result(string found, string degree, string mod, string by, string window, bool table)
        {
            SendStatus();
            string path = Directory.GetCurrentDirectory();

            string[] founds = Directory.GetFiles(path + "\\Base\\" + found);
            string[] degrees = Directory.GetFiles(path + "\\Exponent\\" + degree);
            string[] mods = Directory.GetFiles(path + "\\Modulus\\" + mod);
            var lengts = GenFunctions.ReadString(Window);
            DirectoryInfo di = CD(Choice, found, degree, mod, Window);
            switch (Choice.Split(' ')[0])
            {
                case "Base":
                    switch (Choice.Split(' ')[1])
                    {
                        case "Window":
                            GenWin2(founds, degrees, mods, di.FullName, "Base\\" + found, "Exponent\\" + degree, "Modulus\\" + mod, Window, table);
                            break;
                        case "Exponent":
                            Gen(founds, degrees, mods, di.FullName, "Base\\" + found, "Exponent\\" + degree, "Modulus\\" + mod, Window, table);
                            break;
                        case "Modulus":
                            Gen(founds, mods, degrees, di.FullName, "Base\\" + found, "Modulus\\" + mod, "Exponent\\" + degree, Window, table);
                            break;
                        default: goto smth_go_wrong;
                    }
                    break;

                case "Exponent":
                    switch (Choice.Split(' ')[1])
                    {
                        case "Window":
                            GenWin2(degrees, founds, mods, di.FullName, "Exponent\\" + degree, "Base\\" + found, "Modulus\\" + mod, Window, table);
                            break;
                        case "Base":
                            Gen(degrees, founds, mods, di.FullName, "Exponent\\" + degree, "Base\\" + found, "Modulus\\" + mod, Window, table);
                            break;
                        case "Modulus":
                            Gen(degrees, mods, founds, di.FullName, "Exponent\\" + degree, "Modulus\\" + mod, "Base\\" + found, Window, table);
                            break;
                        default: goto smth_go_wrong;
                    }
                    break;

                case "Modulus":
                    switch (Choice.Split(' ')[1])
                    {
                        case "Window":
                            GenWin2(mods, degrees, founds, di.FullName, "Modulus\\" + mod, "Exponent\\" + degree, "Base\\" + found, Window, table);
                            break;
                        case "Base":
                            Gen(mods, founds, degrees, di.FullName, "Modulus\\" + mod, "Base\\" + found, "Exponent\\" + degree, Window, table);
                            break;
                        case "Exponent":
                            Gen(mods, degrees, founds, di.FullName, "Modulus\\" + mod, "Exponent\\" + degree, "Base\\" + found, Window, table);
                            break;
                        default: goto smth_go_wrong;
                    }

                    break;


                case "Window":
                    switch (Choice.Split(' ')[1])
                    {
                        case "Exponent":
                            GenWin(degrees, founds, mods, di.FullName, "Exponent\\" + degree, "Base\\" + found, "Modulus\\" + mod, Window, table);
                            break;
                        case "Base":
                            GenWin(founds, degrees, mods, di.FullName, "Base\\" + found, "Exponent\\" + degree, "Modulus\\" + mod, Window, table);
                            break;
                        case "Modulus":
                            GenWin(mods, degrees, founds, di.FullName, "Modulus\\" + mod, "Exponent\\" + degree, "Base\\" + found, Window, table);
                            break;
                        default: goto smth_go_wrong;
                    }
                    break;

                default:
                smth_go_wrong:
                    MessageBox.Show("Wrong parametrs", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                    break;
            }
            SendStatus();
        }

        protected double Watch(string Found_file, string Degree_file, string Mod_file, int w, out double table_time)
        {
            Stopwatch stw = new Stopwatch();
            List<double> sums = new List<double>();
            double f_sum = 0;
            double t_d, t_m, t_t;
            table_time = 0;
            int f_count = 0;
            using (StreamReader srf = new StreamReader(Found_file))
            {
                string found;
                while ((found = srf.ReadLine()) != null)
                {
                    double d_sum = 0;
                    int d_count = 0;
                    t_d = 0;
                    using (StreamReader srd = new StreamReader(Degree_file))
                    {
                        string degree;
                        while ((degree = srd.ReadLine()) != null)
                        {
                            double m_sum = 0;
                            int m_count = 0;
                            t_m = 0;
                            using (StreamReader srm = new StreamReader(Mod_file))
                            {
                                string mod;
                                while ((mod = srm.ReadLine()) != null)
                                {
                                    stw.Start();
                                    LoopFunc(BigInteger.Parse(found), BigInteger.Parse(degree), BigInteger.Parse(mod), w, out t_t);
                                    stw.Stop();

                                    m_sum += stw.Elapsed.TotalMilliseconds;
                                    m_count++;
                                    t_m += t_t;
                                    stw.Reset();
                                }

                            }
                            d_sum += m_sum / m_count;
                            t_d += t_m / m_count;
                            d_count++;
                        }
                    }
                    f_sum += d_sum / d_count;
                    table_time += t_d / d_count;
                    f_count++;
                }
            }
            table_time = table_time / f_count;
            return f_sum / f_count;
        }
    }

    #region Binary
    class Bonus1 : BinaryPow
    {
        public Bonus1(string found, string degree, string mod, string by) : base(found, degree, mod, by) { }

        public override string Name() { return "Bonus 1"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod)
        {
             PowFunctions.Bonus1(found, pow, mod); }
    }
    class Bonus2 : BinaryPow
    {
        public Bonus2(string found, string degree, string mod, string by) : base(found, degree, mod, by) { }

        public override string Name() { return "Bonus 2"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod) { PowFunctions.Bonus2(found, pow, mod);  }
    }

    class BinaryRL : BinaryPow
    {
        public BinaryRL(string found, string degree, string mod, string by) : base(found, degree, mod, by){ }
        
        public override string Name() { return "Binary RL"; }

        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod)
        { PowFunctions.BinaryRL(found, pow, mod);  }
    }
    class PowCSharp : BinaryPow
    {
        public PowCSharp(string found, string degree, string mod, string by) : base(found, degree, mod, by) { }

        public override string Name() { return "Pow C#"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod) { BigInteger.ModPow(found, pow, mod);}
    }
    
    class BinaryLR : BinaryPow
    {
        public BinaryLR(string found, string degree, string mod, string by) : base(found, degree, mod, by) { }

        public override string Name() { return "Binary LR"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod) { PowFunctions.BinaryLR(found, pow, mod); }
    }
    class NAFBinaryRL : BinaryPow
    {
        public NAFBinaryRL(string found, string degree, string mod, string by) : base(found, degree, mod, by){ }

        public override string Name() { return "NAF Binary RL"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod) { PowFunctions.NAFBinaryRL(found, pow, mod); }
    }
    class NAFBinaryLR : BinaryPow
    {
        public NAFBinaryLR(string found, string degree, string mod, string by) : base(found, degree, mod, by){ }

        public override string Name() { return "NAF Binary LR"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod) { PowFunctions.NAFBinaryLR(found, pow, mod); }
    }
    class AddSubRL : BinaryPow
    {
        public AddSubRL(string found, string degree, string mod, string by) : base(found, degree, mod, by){ }

        public override string Name() { return "Add Sub RL"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod) { PowFunctions.AddSubRL(found, pow, mod); }
    }
    class AddSubLR : BinaryPow
    {
        public AddSubLR(string found, string degree, string mod, string by) : base(found, degree, mod, by){ }

        public override string Name() { return "Add Sub LR"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod) { PowFunctions.AddSubLR(found, pow, mod); }
    }
    class Joye_double_and_add : BinaryPow
    {
        public Joye_double_and_add(string found, string degree, string mod, string by) : base(found, degree, mod, by){ }

        public override string Name() { return "Joye double and add"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod) { PowFunctions.Joye_double_and_add(found, pow, mod); }
    }
    class MontgomeryLadder : BinaryPow
    {
        public MontgomeryLadder(string found, string degree, string mod, string by) : base(found, degree, mod, by){ }

        public override string Name() { return "Montgomery Ladder"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod) { PowFunctions.MontgomeryLadder(found, pow, mod); }
    }
    class DBNS2RL : BinaryPow
    {
        public DBNS2RL(string found, string degree, string mod, string by) : base(found, degree, mod, by){ }

        public override string Name() { return "DBNS 2 RL"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod) { PowFunctions.DBNS2RL(found, pow, mod); }
    }
    class DBNS2LR : BinaryPow
    {
        public DBNS2LR(string found, string degree, string mod, string by) : base(found, degree, mod, by){ }

        public override string Name() { return "DBNS 2 LR"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod) { PowFunctions.DBNS2LR(found, pow, mod); }
    }
    class DBNS1RL : BinaryPowUp
    {
        bool convert_method;

        public DBNS1RL(string found, string degree, string mod, string by, bool choice, string a_max, string b_max) : base(found, degree, mod, by, a_max, b_max)
        {
            convert_method = choice;
        }

        public override string Name() { return convert_method ? "DBNS 1.1 RL" : "DBNS 1.2 RL"; }

        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, params BigInteger[] list)
        {
            PowFunctions.DBNS1RL(found, pow, mod, convert_method, (int)list[0], (int)list[1]);
        }
    }
    class DBNS1LR : BinaryPowUp
    {
        bool convert_method;

        public DBNS1LR(string found, string degree, string mod, string by, bool choice, string a_max, string b_max) : base(found, degree, mod, by, a_max, b_max)
        {
            convert_method = choice;
        }

        public override string Name() { return convert_method ? "DBNS 1.1 LR" : "DBNS 1.2 LR"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, params BigInteger[] list) { PowFunctions.DBNS1LR(found, pow, mod, convert_method, (int)list[0], (int)list[1]); }
    }
    #endregion
    #region Window
    class WindowRL : WindowPow
    {
        public WindowRL(string found, string degree, string mod, string by, string window, bool table) :
                   base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window RL"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowRL(found, pow, mod, w, out TableTime); }
    }
    class WindowRL_Dic : WindowPow
    {
        public WindowRL_Dic(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window RL_Dic"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowRL_Dic(found, pow, mod, w, out TableTime); }
    }
    class WindowLR : WindowPow
    {
        public WindowLR(string found, string degree, string mod, string by, string window, bool table) :
                   base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LR"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLR(found, pow, mod, w, out TableTime); }
    }
    class WindowLR_Dic : WindowPow
    {
        public WindowLR_Dic(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LR_Dic"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLR_Dic(found, pow, mod, w, out TableTime); }
    }
    class SlidingRL : WindowPow
    {
        public SlidingRL(string found, string degree, string mod, string by, string window, bool table) :
                   base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Sliding Window RL"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.SlidingRL(found, pow, mod, w, out TableTime); }
    }
    class SlidingRL_Dic : WindowPow
    {
        public SlidingRL_Dic(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Sliding Window RL_Dic"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.SlidingRL_Dic(found, pow, mod, w, out TableTime); }
    }
    class SlidingLR : WindowPow
    {
        public SlidingLR(string found, string degree, string mod, string by, string window, bool table) :
                   base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Sliding Window LR"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.SlidingLR(found, pow, mod, w, out TableTime); }
    }
    class SlidingLR_Dic : WindowPow
    {
        public SlidingLR_Dic(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Sliding Window LR_Dic"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.SlidingLR_Dic(found, pow, mod, w, out TableTime); }
    }
    class NAFSlidingRL : WindowPow
    {
        public NAFSlidingRL(string found, string degree, string mod, string by, string window, bool table) :
                   base(found, degree, mod, by, window, table) { }

        public override string Name() { return "NAF Sliding RL"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.NAFSlidingRL(found, pow, mod, w, out TableTime); }
    }
    class NAFSlidingLR : WindowPow
    {
        public NAFSlidingLR(string found, string degree, string mod, string by, string window, bool table) :
                   base(found, degree, mod, by, window, table) { }

        public override string Name() { return "NAF Sliding LR"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.NAFSlidingLR(found, pow, mod, w, out TableTime); }
    }
    class NAFWindowRL : WindowPow
    {
        public NAFWindowRL(string found, string degree, string mod, string by, string window, bool table) :
                   base(found, degree, mod, by, window, table) { }

        public override string Name() { return "NAF Window RL"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.NAFWindowRL(found, pow, mod, w, out TableTime); }
    }
    class NAFWindowLR : WindowPow
    {
        public NAFWindowLR(string found, string degree, string mod, string by, string window, bool table) :
                   base(found, degree, mod, by, window, table) { }

        public override string Name() { return "NAF Window LR"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.NAFWindowLR(found, pow, mod, w, out TableTime); }
    }
    class wNAFSlidingRL : WindowPow
    {
        public wNAFSlidingRL(string found, string degree, string mod, string by, string window, bool table) :
                   base(found, degree, mod, by, window, table) { }

        public override string Name() { return "wNAF Sliding RL"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.wNAFSlidingRL(found, pow, mod, w, out TableTime); }
    }
    class wNAFSlidingLR : WindowPow
    {
        public wNAFSlidingLR(string found, string degree, string mod, string by, string window, bool table) :
                   base(found, degree, mod, by, window, table) { }

        public override string Name() { return "wNAF Sliding LR"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.wNAFSlidingLR(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod1 : WindowPow
    {
        public WindowLRMod1(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod1"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod1(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod2 : WindowPow
    {
        public WindowLRMod2(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod2"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod2(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod3 : WindowPow
    {
        public WindowLRMod3(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod3"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod3(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod : WindowPow
    {
        public WindowLRMod(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod1_Shift : WindowPow
    {
        public WindowLRMod1_Shift(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod1_Shift"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod1_Shift(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod2_Shift : WindowPow
    {
        public WindowLRMod2_Shift(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod2_Shift"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod2_Shift(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod3_Shift : WindowPow
    {
        public WindowLRMod3_Shift(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod3_Shift"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod3_Shift(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod_Shift : WindowPow
    {
        public WindowLRMod_Shift(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod_Shift"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod_Shift(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod1_Upgrade : WindowPow
    {
        public WindowLRMod1_Upgrade(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod1_Upgrade"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod1_Upgrade(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod2_Upgrade : WindowPow
    {
        public WindowLRMod2_Upgrade(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod2_Upgrade"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod2_Upgrade(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod3_Upgrade : WindowPow
    {
        public WindowLRMod3_Upgrade(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod3_Upgrade"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod3_Upgrade(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod_Upgrade : WindowPow
    {
        public WindowLRMod_Upgrade(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod_Upgrade"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod_Upgrade(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod1_NoBinary : WindowPow
    {
        public WindowLRMod1_NoBinary(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod1_NoBinary"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod1_Upgrade(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod2_NoBinary : WindowPow
    {
        public WindowLRMod2_NoBinary(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod2_NoBinary"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod2_NoBinary(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod3_NoBinary : WindowPow
    {
        public WindowLRMod3_NoBinary(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod3_NoBinary"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod3_NoBinary(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod_NoBinary : WindowPow
    {
        public WindowLRMod_NoBinary(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod_NoBinary"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod_NoBinary(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod1_Final : WindowPow
    {
        public WindowLRMod1_Final(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod1_Final"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod1_Upgrade(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod2_Final : WindowPow
    {
        public WindowLRMod2_Final(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod2_Final"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod2_Final(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod3_Final : WindowPow
    {
        public WindowLRMod3_Final(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod3_Final"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod3_Final(found, pow, mod, w, out TableTime); }
    }
    class WindowLRMod_Final : WindowPow
    {
        public WindowLRMod_Final(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Window LRMod_Final"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.WindowLRMod_Final(found, pow, mod, w, out TableTime); }
    }
    class Sliding_Prime : WindowPow
    {
        private bool state;
        public Sliding_Prime(string found, string degree, string mod, string by, string window, bool table, bool state) :
            base(found, degree, mod, by, window, table)
        {
            this.state = state;
        }

        public override string Name() { return state ? "Sliding Prime with_11_11" : "Sliding Prime"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.Sliding_Prime(found, pow, mod, w, out TableTime, state); }
    }
    class Adaptive : WindowPow
    {
        public Adaptive(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Adaptive"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.Adaptive(found, pow, mod, w, out TableTime); }
    }
    class BonusWindow : WindowPow
    {
        public BonusWindow(string found, string degree, string mod, string by, string window, bool table) :
            base(found, degree, mod, by, window, table) { }

        public override string Name() { return "Bonus Window"; }
        protected override void LoopFunc(BigInteger found, BigInteger pow, BigInteger mod, int w, out double TableTime) { BigInteger e = PowFunctions.Bonus(found, pow, mod, w, out TableTime); }
    }
    #endregion
}
