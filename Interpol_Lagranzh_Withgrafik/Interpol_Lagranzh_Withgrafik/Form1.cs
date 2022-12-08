using System;
using static System.Math;
using static System.Console;
using static System.Convert;
using System.Windows.Forms;
using ZedGraph;
using System.Drawing;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using TestMySpline;
using ZedGraph;
namespace Lagrange
{
    public partial class Form1 : Form
    {
        private static Func<double, double> f;
        private static Func<double, double> dfx;
        private static Func<double, double> f1 = x => Abs(x);
        private static Func<double, double> f2 = x => x*x-2*x+1;
        private static Func<double, double> f3 = x => Exp(-1 * Pow(x, 2));
        private static Func<double, double> df1 = x => Abs(x)/x;
        private static Func<double, double> df2 = x => 2*x-2;
        private static Func<double, double>df3 = x =>(-2*x)/ Exp(Pow(x, 2));
        public Form1()
        {
            InitializeComponent();
            DrawGraph();
        }


        private static double Lagrange(double a, ref double[] x, ref double[] y)
        {
            double f = 0;
            for (int i = 0; i < x.Length; ++i)
            {
                double l = 1;
                for (int j = 0; j < x.Length; ++j)
                    if (i != j)
                        l *= (a - x[j]) / (x[i] - x[j]);
                l *= y[i];
                f += l;
            }
            return f;
        }




        SplineTuple[] splines; // Сплайн

        // Структура, описывающая сплайн на каждом сегменте сетки
        private struct SplineTuple
        {
            public double a, b, c, d, x;
        }

        // Построение сплайна
        // x - узлы сетки, должны быть упорядочены по возрастанию, кратные узлы запрещены
        // y - значения функции в узлах сетки
        // n - количество узлов сетки
        public void BuildSpline(double[] x, double[] y, int n)
        {
            // Инициализация массива сплайнов
            splines = new SplineTuple[n];
            for (int i = 0; i < n; ++i)
            {
                splines[i].x = x[i];
                splines[i].a = y[i];
            }
            splines[0].c = splines[n - 1].c = 0.0;

            // Решение СЛАУ относительно коэффициентов сплайнов c[i] методом прогонки для трехдиагональных матриц
            // Вычисление прогоночных коэффициентов - прямой ход метода прогонки
            double[] alpha = new double[n - 1];
            double[] beta = new double[n - 1];
            alpha[0] = beta[0] = 0.0;
            for (int i = 1; i < n - 1; ++i)
            {
                double hi = x[i] - x[i - 1];
                double hi1 = x[i + 1] - x[i];
                double A = hi;
                double C = 2.0 * (hi + hi1);
                double B = hi1;
                double F = 6.0 * ((y[i + 1] - y[i]) / hi1 - (y[i] - y[i - 1]) / hi);
                double z = (A * alpha[i - 1] + C);
                alpha[i] = -B / z;
                beta[i] = (F - A * beta[i - 1]) / z;
            }

            // Нахождение решения - обратный ход метода прогонки
            for (int i = n - 2; i > 0; --i)
            {
                splines[i].c = alpha[i] * splines[i + 1].c + beta[i];
            }

            // По известным коэффициентам c[i] находим значения b[i] и d[i]
            for (int i = n - 1; i > 0; --i)
            {
                double hi = x[i] - x[i - 1];
                splines[i].d = (splines[i].c - splines[i - 1].c) / hi;
                splines[i].b = hi * (2.0 * splines[i].c + splines[i - 1].c) / 6.0 + (y[i] - y[i - 1]) / hi;
            }
        }

        // Вычисление значения интерполированной функции в произвольной точке
        public double Interpolate(double x)
        {
            if (splines == null)
            {
                return double.NaN; // Если сплайны ещё не построены - возвращаем NaN
            }

            int n = splines.Length;
            SplineTuple s;

            if (x <= splines[0].x) // Если x меньше точки сетки x[0] - пользуемся первым эл-тов массива
            {
                s = splines[0];
            }
            else if (x >= splines[n - 1].x) // Если x больше точки сетки x[n - 1] - пользуемся последним эл-том массива
            {
                s = splines[n - 1];
            }
            else // Иначе x лежит между граничными точками сетки - производим бинарный поиск нужного эл-та массива
            {
                int i = 0;
                int j = n - 1;
                while (i + 1 < j)
                {
                    int k = i + (j - i) / 2;
                    if (x <= splines[k].x)
                    {
                        j = k;
                    }
                    else
                    {
                        i = k;
                    }
                }
                s = splines[j];
            }

            double dx = x - s.x;
            // Вычисляем значение сплайна в заданной точке по схеме Горнера (в принципе, "умный" компилятор применил бы схему Горнера сам, но ведь не все так умны, как кажутся)
            return s.a + (s.b + (s.c / 2.0 + s.d * dx / 6.0) * dx) * dx;
        }
        //**************************************************************************************************************



        static public double Newton(double x, int n, double[] MasX, double[] MasY, double step)
        {
            double[,] mas = new double[n + 2, n + 1];
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < n + 1; j++)
                {
                    if (i == 0)
                        mas[i, j] = MasX[j];
                    else if (i == 1)
                        mas[i, j] = MasY[j];
                }
            }
            int m = n;
            for (int i = 2; i < n + 2; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    mas[i, j] = mas[i - 1, j + 1] - mas[i - 1, j];
                }
                m--;
            }

            double[] dy0 = new double[n + 1];

            for (int i = 0; i < n + 1; i++)
            {
                dy0[i] = mas[i + 1, 0];
            }

            double res = dy0[0];
            double[] xn = new double[n];
            xn[0] = x - mas[0, 0];

            for (int i = 1; i < n; i++)
            {
                double ans = xn[i - 1] * (x - mas[0, i]);
                xn[i] = ans;
                ans = 0;
            }

            int m1 = n + 1;
            int fact = 1;
            for (int i = 1; i < m1; i++)
            {
                fact = fact * i;
                res = res + (dy0[i] * xn[i - 1]) / (fact * Math.Pow(step, i));
            }

            return res;
        }

        public void drownio() // метод рисования 
        {


            if (radioButton1.Checked)
            {
                f = f1; dfx = df1;
            }
            else if (radioButton2.Checked) { f = f2; dfx = df2; }
            else if (radioButton3.Checked) { f = f3; dfx = df3; }
            double[] x = new double[1000];
            double[] y = new double[1000];
            double Xmin = Convert.ToDouble(textBox2.Text);
            double Xmax = Convert.ToDouble(textBox3.Text);
            int X = Convert.ToInt32(textBox1.Text);
            double  h = (Xmax - Xmin) / (X - 1);
            int k = X;
            for (int i = 0; i < k; i++)
            {
                x[i] = Xmin + i * h;
                y[i] = f(x[i]);
            }
           
            GraphPane grap = zedGraph.GraphPane; // графический объект
            grap.CurveList.Clear(); // очищаю график
            
            PointPairList point = new PointPairList();
            // шаг для равноотстоящих узлов
            for (double i = Xmin; i <= Xmax; i += h)
            {
                point.Add(i, f(i)); // заполняю точки
            }

            var pane = zedGraph.MasterPane;
            PointPairList point2 = new PointPairList();
            for (double i = Xmin; i <=Xmax; i += h)
            {
                point2.Add(i,Newton(i,X,x,y,h));
            }



            //LineItem line =  grap.AddCurve("График", point, Color.Red, SymbolType.None);
            //grap.AddCurve("Ньютона", point2, Color.Blue, SymbolType.None);
            LineItem line = grap.AddCurve("Ньютона", point2, Color.Blue, SymbolType.None);
            zedGraph.AxisChange();
            zedGraph.Invalidate();
        }










        //**************************************************************************************************************

        public void drow() // метод рисования 
        {
            if (radioButton1.Checked) f = f1;
            else if (radioButton2.Checked) f = f2;
            else if (radioButton3.Checked) f = f3;
            double X = Convert.ToDouble(textBox1.Text);
            GraphPane grap = zedGraph.GraphPane; // графический объект
            grap.CurveList.Clear(); // очищаю график
            double Xmin = Convert.ToDouble(textBox2.Text);
            double Xmax = Convert.ToDouble(textBox3.Text);
            PointPairList point = new PointPairList();
            double h = 1.0 / X;
            for (double i = Xmin; i < Xmax; i += h)
            {
                point.Add(i, f(i)); // заполняю точки
            }

            var pane = zedGraph.MasterPane;
            var toInterpolate = new List<float>();
            for (double i = Xmin; i < Xmax; i += 0.01)
            {
                toInterpolate.Add((float)i);
            }
            var spline = new CubicSpline();
            var ys = spline.FitAndEval(point.Select(p => (float)p.X).ToArray(),
                point.Select(p => (float)p.Y).ToArray(),
                toInterpolate.ToArray());

            var splinePairs = new PointPairList(toInterpolate.Select(e => (double)e).ToArray(), ys.Select(e => (double)e).ToArray());

            LineItem line = grap.AddCurve("График", point, Color.Red, SymbolType.None); // стою линию 
            grap.AddCurve("Сплайнами", splinePairs, Color.Blue, SymbolType.None);
            zedGraph.AxisChange();
            zedGraph.Invalidate();
        }

        private void DrawGraph()
        {
            WriteLine("Задайте границы a и b: ");
            double a, b, h;
            a = Convert.ToInt32(textBox2.Text);
            b = Convert.ToInt32(textBox3.Text);
            int n, m;
            WriteLine("Задайте кол-во точек n: ");
            n = Convert.ToInt32(textBox1.Text);
            h = (b - a) / (n - 1); // шаг для равноотстоящих узлов
            if (radioButton1.Checked) f = f1;
            else if (radioButton2.Checked) f = f2;
            else if (radioButton3.Checked) f = f3;
            double[] x = new double[n];
            double[] y = new double[n];
            for (int i = 0; i < n; i++)
            {
                if (radioButton4.Checked) { x[i] = a + i * h; } // равноотстоящие узлы
                if (radioButton5.Checked) { x[i] = Cos((2 * i + 1) * PI / (2 * n)); } // узлы Чебышева
                y[i] = f(x[i]);
            }
            int k;
            WriteLine("Задайте кол-во точек k: ");
            //k = ToInt32(ReadLine());
            k = 1000;
            double h1, h2;
            h1 = h2 = (b - a) / (k - 1);
            double[] x1 = new double[k];
            double[] y1 = new double[k];
            for (int i = 0; i < k; i++)
            {
                x1[i] = a + i * h1;
                y1[i] = Lagrange(x1[i], ref x, ref y);
            }
            double[] x2 = new double[k];
            double[] y2 = new double[k];
            for (int i = 0; i < k; i++)
            {
                x2[i] = a + i * h2;
                y2[i] = f(x2[i]);
            }

            // Получим панель для рисования
            GraphPane pane = zedGraph.GraphPane;

            // Очистим список кривых на тот случай, если до этого сигналы уже были нарисованы
            pane.CurveList.Clear();

            // Создадим список точек для кривой f1(x)
            PointPairList f1_list = new PointPairList();

            // Создадим список точек для кривой f2(x)
            PointPairList f2_list = new PointPairList();

            // !!!
            // Заполним массив точек для кривой f1(x)
            for (int i = 0; i < x1.Length; ++i)
            {
                f1_list.Add(x1[i], y1[i]);
            }

            // !!!
            // Заполним массив точек для кривой f2(x)
            // Интервал и шаги по X могут не совпадать на разных кривых
            for (int i = 0; i < x1.Length; ++i)
            {
                f2_list.Add(x2[i], y2[i]);
            }

            // !!!
            // Создадим кривую с названием "Lagrange",
            // которая будет рисоваться голубым цветом (Color.Blue),
            // Опорные точки выделяться не будут (SymbolType.None)
            pane.AddCurve("Lagrange", f1_list, Color.Blue, SymbolType.None);

            // !!!
            // Создадим кривую с названием "Function",
            // которая будет рисоваться красным цветом (Color.Red),
            // Опорные точки выделяться не будут (SymbolType.None)
            pane.AddCurve("Function", f2_list, Color.Red, SymbolType.None);

            // Вызываем метод AxisChange (), чтобы обновить данные об осях.
            // В противном случае на рисунке будет показана только часть графика,
            // которая умещается в интервалы по осям, установленные по умолчанию
            zedGraph.AxisChange();

            // Изменим текст заголовка графика
            pane.Title.Text = "Графики многочлена и функции";
            // Изменим тест надписи по оси x
            pane.XAxis.Title.Text = "x";
            // Изменим тест надписи по оси y
            pane.YAxis.Title.Text = "y";

            // Обновляем график
            zedGraph.Invalidate();
        }
        private void button1_Click(object sender, EventArgs e)
        {
            if (radioButton7.Checked) { DrawGraph(); }
            if (radioButton6.Checked) { drow(); }
            if (radioButton8.Checked) {
                radioButton4.Checked = true;
                drownio(); }

        }

        private void label2_Click(object sender, EventArgs e)
        {

        }
    }
}
