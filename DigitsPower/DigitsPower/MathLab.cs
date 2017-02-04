//using System;
//using System.Collections.Generic;
//using System.Windows.Forms;
//using System.Numerics;
//using System.IO;
//using System.Diagnostics;
//using System.Threading.Tasks;
//using System.Threading;

////using System.Security.Cryptography;

//namespace DigitsPower
//{
//    //1.Binary RL
//    //2.Binary LR
//    //3.Window RL
//    //4.Window LR
//    //5.Slide Window RL
//    //6.Slide Window LR
//    //7.1.NAF Binary RL
//    //7.2.NAF Binary RL_2
//    //8.NAF Binary LR
//    //9.NAF Slide RL
//    //10.NAF Slide LR
//    //11.NAF Window RL
//    //12.NAF Window LR
//    //13.wNAF Slide RL
//    //14.wNAF Slide LR
//    //15.Add Sub RL
//    //16.Add Sub LR
//    //17.Joye double & add
//    //18.Montgomery ladder
//    //19.DBNS 1 RL
//    //20.DBNS 1 LR
//    //21.DBNS 2 RL
//    //22.DBNS 2 LR

//    public partial class MainForm : Form
//    {
//        string choice;
//        string choicew;
//        public MainForm()
//        {
//            InitializeComponent();
//            CreateDirectories();
//            Axis1Box.SelectedIndex = 0;
//            Axis2Box.SelectedIndex = 1;
//            tabControl1.SelectedIndex = 2;
//        }

//        #region InitializeMethods
//        private void CreateDirectories()
//        {
//            string path = Directory.GetCurrentDirectory();

//            try
//            {
//                Directory.CreateDirectory(path + "\\Foundations");
//                Directory.CreateDirectory(path + "\\Degrees");
//                Directory.CreateDirectory(path + "\\Mods");
//                Directory.CreateDirectory(path + "\\Results");
//                UpdateDirectoryList();
//                OperCheckList.SetItemChecked(0, true);
//                OperCheckList.SetItemChecked(1, true);
//                if (FoundDir.Items.Count != 0)
//                {
//                    FoundDir.SelectedIndex = 0;
//                }
//                if (DegreeDir.Items.Count != 0)
//                {
//                    DegreeDir.SelectedIndex = 0;
//                }
//                if (ModsDir.Items.Count != 0)
//                {
//                    ModsDir.SelectedIndex = 0;
//                }

//            }
//            catch (Exception err)
//            {
//                MessageBox.Show(err.Message, "Warning", MessageBoxButtons.OK, MessageBoxIcon.Warning);
//            }
//        }

//        private void UpdateDirectoryList()
//        {
//            try
//            {
//                string path = Directory.GetCurrentDirectory();
//                string[] dirs = Directory.GetDirectories(path + "\\Degrees");
//                DegreesList.Items.Clear();
//                DegreeDir.Items.Clear();
//                foreach (string s in dirs)
//                {
//                    string[] dirs2 = s.Split('\\');
//                    DegreesList.Items.Add(dirs2[dirs2.Length - 1]);
//                    DegreeDir.Items.Add(dirs2[dirs2.Length - 1]);
//                }
//                dirs = Directory.GetDirectories(path + "\\Foundations");
//                FoundationsList.Items.Clear();
//                FoundDir.Items.Clear();
//                foreach (string s in dirs)
//                {
//                    string[] dirs2 = s.Split('\\');
//                    FoundationsList.Items.Add(dirs2[dirs2.Length - 1]);
//                    FoundDir.Items.Add(dirs2[dirs2.Length - 1]);
//                }

//                dirs = Directory.GetDirectories(path + "\\Mods");
//                ModsList.Items.Clear();
//                ModsDir.Items.Clear();
//                foreach (string s in dirs)
//                {
//                    string[] dirs2 = s.Split('\\');
//                    ModsList.Items.Add(dirs2[dirs2.Length - 1]);
//                    ModsDir.Items.Add(dirs2[dirs2.Length - 1]);
//                }
//                dirs = Directory.GetDirectories(path + "\\Results");
//            }
//            catch (Exception err)
//            {
//                MessageBox.Show(err.Message, "Warning", MessageBoxButtons.OK, MessageBoxIcon.Warning);
//            }
//        }
//        #endregion

//        #region GenerationMethods
//        private void GenFile(GenFunctions.random_num func, string dir, string len, string count, string type, string radix)
//        {
//            string path = Directory.GetCurrentDirectory();
//            var lenght = GenFunctions.ReadString(len);
//            try
//            {
//                var di = Directory.CreateDirectory($"{path}\\{dir}\\{type}{lenght[0]}_{lenght[lenght.Count - 1]}_{count}({radix})#{DateTime.Now.ToLocalTime().ToString().Replace(':', '-')}");
//                FileStream fin;
//                for (int j = 0; j < lenght.Count; j++)
//                {
//                    path = di.FullName + "\\" + lenght[j] + ".txt";

//                    fin = new FileStream(path, FileMode.Create, FileAccess.Write);
//                    // Open the stream and read it back.
//                    using (StreamWriter sr = new StreamWriter(fin))
//                    {
//                        for (int i = 0; i < Int32.Parse(count); i++)
//                        {
//                            sr.WriteLine(func(lenght[j], radix));
//                        }
//                    }
//                    fin.Close();
//                }
//            }
//            catch (Exception err)
//            {
//                MessageBox.Show(err.Message, "Warning", MessageBoxButtons.OK, MessageBoxIcon.Warning);
//            }
//            UpdateDirectoryList();
//        }

//        private void GenFound_Click(object sender, EventArgs e)
//        {
//            GenFile(GenFunctions.random_max, "Foundations", FoundLenght.Text, FoundCount.Text, "Found", comboBox2.Text);
//        }

//        private void GenDegree_Click(object sender, EventArgs e)
//        {
//            GenFile(GenFunctions.random_max, "Degrees", DegreeLenght.Text, DegreeCount.Text, "Degree", comboBox3.Text);
//        }

//        private void GenMod_Click(object sender, EventArgs e)
//        {
//            switch (ModType.Text)
//            {
//                case "Degree of 2":
//                    GenFile(GenFunctions.random_two, "Mods", ModLenght.Text, ModCount.Text, "Mod_Degree of two_", comboBox4.Text);
//                    break;
//                case "Odd number":
//                    GenFile(GenFunctions.random_odd, "Mods", ModLenght.Text, ModCount.Text, "Mod_Odd number_", comboBox4.Text);
//                    break;
//                case "Prime number":
//                    GenFile(GenFunctions.random_simple, "Mods", ModLenght.Text, ModCount.Text, "Mod_Prime number_", comboBox4.Text);
//                    break;
//                default:
//                    GenFile(GenFunctions.random_max, "Mods", ModLenght.Text, ModCount.Text, "Mod_Random number_", comboBox4.Text);
//                    break;
//            }

//        }
//        #endregion

//        #region CalculateMethods
//        private void ResultsButton_Click(object sender, EventArgs e)
//        {
//            ResultsButton.Enabled = false;
//            AdditionalParameters.A = long.Parse(textBox5.Text);
//            AdditionalParameters.B = long.Parse(textBox4.Text);

//            choice = comboBox1.Text;

//            if (Axis1Box.SelectedItem == Axis2Box.SelectedItem)
//            {
//                MessageBox.Show("You must select two different axis!", "Warning", MessageBoxButtons.OK, MessageBoxIcon.Warning);
//                return;
//            }
//            choicew = Axis1Box.SelectedItem + " " + Axis2Box.SelectedItem;
//            if (OperCheckList.CheckedIndices.Count > 0)
//            {
//                (new Thread(GetResults)).Start();
//            }
//        }

//        private void GetResults()
//        {
//            for (int i = 0; i < OperCheckList.CheckedIndices.Count; i++)
//            {
//                #region Binary
//                if (OperCheckList.CheckedIndices[i] == 0) { (new BinaryRL(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choice)).Create_Result(); }
//                if (OperCheckList.CheckedIndices[i] == 1) { (new BinaryLR(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choice)).Create_Result(); }
//                if (OperCheckList.CheckedIndices[i] == 6) { (new NAFBinaryRL(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choice)).Create_Result(); }
//                if (OperCheckList.CheckedIndices[i] == 7) { (new NAFBinaryLR(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choice)).Create_Result(); }
//                if (OperCheckList.CheckedIndices[i] == 14) { (new AddSubRL(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choice)).Create_Result(); }
//                if (OperCheckList.CheckedIndices[i] == 15) { (new AddSubLR(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choice)).Create_Result(); }
//                if (OperCheckList.CheckedIndices[i] == 16) { (new Joye_double_and_add(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choice)).Create_Result(); }
//                if (OperCheckList.CheckedIndices[i] == 17) { (new MontgomeryLadder(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choice)).Create_Result(); }
//                if (OperCheckList.CheckedIndices[i] == 18) { (new DBNS1RL(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choice, true)).Create_Result(); }
//                if (OperCheckList.CheckedIndices[i] == 19) { (new DBNS1RL(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choice, false)).Create_Result(); }
//                if (OperCheckList.CheckedIndices[i] == 20) { (new DBNS1LR(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choice, true)).Create_Result(); }
//                if (OperCheckList.CheckedIndices[i] == 21) { (new DBNS1LR(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choice, false)).Create_Result(); }
//                if (OperCheckList.CheckedIndices[i] == 22) { (new DBNS2RL(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choice)).Create_Result(); }
//                if (OperCheckList.CheckedIndices[i] == 23) { (new DBNS2RL(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choice)).Create_Result(); }
//                #endregion

//                #region Window
//                if (OperCheckList.CheckedIndices[i] == 2) { (new WindowRL()).Create_Result(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choicew, WinMode.Text, TableWith.Checked); }
//                if (OperCheckList.CheckedIndices[i] == 3) { (new WindowLR()).Create_Result(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choicew, WinMode.Text, TableWith.Checked); }
//                if (OperCheckList.CheckedIndices[i] == 4) { (new SlideRL()).Create_Result(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choicew, WinMode.Text, TableWith.Checked); }
//                if (OperCheckList.CheckedIndices[i] == 5) { (new SlideRL()).Create_Result(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choicew, WinMode.Text, TableWith.Checked); }
//                if (OperCheckList.CheckedIndices[i] == 8) { (new NAFSlideRL()).Create_Result(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choicew, WinMode.Text, TableWith.Checked); }
//                if (OperCheckList.CheckedIndices[i] == 9) { (new NAFSlideLR()).Create_Result(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choicew, WinMode.Text, TableWith.Checked); }
//                if (OperCheckList.CheckedIndices[i] == 10) { (new NAFWindowRL()).Create_Result(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choicew, WinMode.Text, TableWith.Checked); }
//                if (OperCheckList.CheckedIndices[i] == 11) { (new NAFWindowLR()).Create_Result(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choicew, WinMode.Text, TableWith.Checked); }
//                if (OperCheckList.CheckedIndices[i] == 12) { (new wNAFSlideRL()).Create_Result(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choicew, WinMode.Text, TableWith.Checked); }
//                if (OperCheckList.CheckedIndices[i] == 13) { (new wNAFSlideLR()).Create_Result(FoundDir.SelectedItem.ToString(), DegreeDir.SelectedItem.ToString(), ModsDir.SelectedItem.ToString(), choicew, WinMode.Text, TableWith.Checked); }
//                #endregion

//            }

//            string path = Directory.GetCurrentDirectory();
//            string[] dirs = Directory.GetDirectories(path + "\\Results");
//            var result = MessageBox.Show("All done", "Result",
//                                         MessageBoxButtons.OK,
//                                         MessageBoxIcon.Information);
//            ResultsButton.Enabled = true;
//        }

//        private void Calculatebutton_Click(object sender, EventArgs e)
//        {
//            AdditionalParameters.A = Int64.Parse(textA.Text);
//            AdditionalParameters.B = Int64.Parse(textB.Text);

//            BigInteger mod = BigInteger.Parse(modText.Text);
//            BigInteger pow = BigInteger.Parse(PowerText.Text);
//            int window = Int32.Parse(WindowText.Text);
//            BigInteger num = BigInteger.Parse(NumberText.Text);
//            List<string> s = new List<string>();
//            double table = 0;
//            OperationsResult.Items.Clear();

//            if (OperationsList.CheckedIndices.Count > 0)
//            {
//                for (int i = 0; i < OperationsList.CheckedIndices.Count; i++)
//                {

//                    if (OperationsList.CheckedIndices[i] == 0) { OperationsResult.Items.Add("Binary RL\t: " + (PowFunctions.BinaryRL(num, pow, mod)).ToString()); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 1) { OperationsResult.Items.Add("Binary LR\t: " + (PowFunctions.BinaryLR(num, pow, mod)).ToString()); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 2) { OperationsResult.Items.Add("Window RL\t: " + (PowFunctions.WindowRL(num, pow, mod, window, out table)).ToString()); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 3) { OperationsResult.Items.Add("Window LR\t: " + (PowFunctions.WindowLR(num, pow, mod, window, out table)).ToString()); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 4) { OperationsResult.Items.Add("Slide RL\t\t: " + (PowFunctions.SlideRL(num, pow, mod, window, out table)).ToString()); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 5) { OperationsResult.Items.Add("Slide LR\t\t: " + (PowFunctions.SlideLR(num, pow, mod, window, out table)).ToString()); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 6) { OperationsResult.Items.Add("NAF Binary RL\t: " + (PowFunctions.NAFBinaryRL(num, pow, mod).ToString())); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 7) { OperationsResult.Items.Add("NAF Binary LR\t: " + (PowFunctions.NAFBinaryLR(num, pow, mod).ToString())); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 8) { OperationsResult.Items.Add("NAF Slide RL\t: " + (PowFunctions.NAFSlideRL(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 9) { OperationsResult.Items.Add("NAF Slide LR\t: " + (PowFunctions.NAFSlideLR(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 10) { OperationsResult.Items.Add("NAF Window RL\t: " + (PowFunctions.NAFWindowRL(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 11) { OperationsResult.Items.Add("NAF Window LR\t: " + (PowFunctions.NAFWindowLR(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 12) { OperationsResult.Items.Add("wNAF Slide RL\t: " + (PowFunctions.wNAFSlideRL(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 13) { OperationsResult.Items.Add("wNAF Slide RL\t: " + (PowFunctions.wNAFSlideLR(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 14) { OperationsResult.Items.Add("Add Sub RL\t: " + (PowFunctions.AddSubRL(num, pow, mod).ToString())); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 15) { OperationsResult.Items.Add("Add Sub LR\t: " + (PowFunctions.AddSubLR(num, pow, mod).ToString())); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 16) { OperationsResult.Items.Add("Joye_double_and_add\t: " + (PowFunctions.Joye_double_and_add(num, pow, mod).ToString())); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 17) { OperationsResult.Items.Add("MontgomeryLadder\t: " + (PowFunctions.MontgomeryLadder(num, pow, mod).ToString())); OperationsResult.Update(); }


//                    if (OperationsList.CheckedIndices[i] == 18) { OperationsResult.Items.Add("DBNS1RL 1\t: " + (PowFunctions.DBNS1RL(num, pow, mod, true).ToString())); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 19) { OperationsResult.Items.Add("DBNS1RL 2\t: " + (PowFunctions.DBNS1RL(num, pow, mod, false).ToString())); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 20) { OperationsResult.Items.Add("DBNS1LR 1\t: " + (PowFunctions.DBNS1LR(num, pow, mod, true).ToString())); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 21) { OperationsResult.Items.Add("DBNS1LR 2\t: " + (PowFunctions.DBNS1LR(num, pow, mod, false).ToString())); OperationsResult.Update(); }

//                    if (OperationsList.CheckedIndices[i] == 22) { OperationsResult.Items.Add("DBNS2RL\t: " + (PowFunctions.DBNS2RL(num, pow, mod, PowFunctions.Euclid_2_1).ToString())); OperationsResult.Update(); }
//                    if (OperationsList.CheckedIndices[i] == 23) { OperationsResult.Items.Add("DBNS2LR\t: " + (PowFunctions.DBNS2LR(num, pow, mod, PowFunctions.Euclid_2_1).ToString())); OperationsResult.Update(); }
//                }
//            }
//        }
//        #endregion

//        #region CheckButtons
//        private void CheckAllbutton_Click(object sender, EventArgs e)
//        {
//            for (int i = 0; i < OperationsList.Items.Count; i++)
//            {
//                OperationsList.SetItemChecked(i, true);
//            }
//        }

//        private void CheckNone_Click(object sender, EventArgs e)
//        {
//            for (int i = 0; i < OperCheckList.Items.Count; i++)
//            {
//                OperCheckList.SetItemChecked(i, false);
//            }
//        }

//        private void CheckNonebutton_Click(object sender, EventArgs e)
//        {
//            for (int i = 0; i < OperationsList.Items.Count; i++)
//            {
//                OperationsList.SetItemChecked(i, false);

//            }
//        }

//        private void button2_Click(object sender, EventArgs e)
//        {
//            for (int i = 0; i < OperCheckList.Items.Count; i++)
//            {
//                OperCheckList.SetItemChecked(i, true);
//            }
//        }

//        #endregion

//        #region OpenDirectories
//        private void ModsDir_DoubleClick(object sender, EventArgs e)
//        {
//            Process.Start(@"~\bin\Debug\Mods");
//        }

//        private void DegreeDir_DoubleClick(object sender, EventArgs e)
//        {
//            Process.Start(@"~\bin\Debug\Degrees");
//        }

//        private void FoundDir_DoubleClick(object sender, EventArgs e)
//        {
//            Process.Start(@"~\bin\Debug\Foudations");
//        }


//        //private void ResultList_DoubleClick(object sender, EventArgs e)
//        //{
//        //    var x = Directory.GetCurrentDirectory();
//        //    var p = x + @"\bin\Debug\Results";
//        //    System.Diagnostics.Process.Start(new System.Diagnostics.ProcessStartInfo()
//        //    {
//        //        FileName = p,
//        //        UseShellExecute = true,
//        //        Verb = "open"
//        //    });
//        //}
//        #endregion
//    }
//    public static class AdditionalParameters
//    {
//        public static long A = 15;
//        public static long B = 17;
//        public static Inverse inv = PowFunctions.Euclid_2_1;

//    }
//}
