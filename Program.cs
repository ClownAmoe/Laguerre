using System;
using System.IO;

public class LaguerreFunc
{
    private int _n;
    private double _t;
    private double _sigma;
    private double _beta;

    public LaguerreFunc(int n, double t, double sigma = 4, double beta = 2)
    {
        N = n;
        T = t;
        Sigma = sigma;
        Beta = beta;
    }

    public int N
    {
        get { return _n; }
        set
        {
            if (value < 0)
            {
                throw new ArgumentException("N cannot be negative.");
            }
            _n = value;
        }
    }

    public double T
    {
        get { return _t; }
        set
        {
            if (value < 0)
            {
                throw new ArgumentException("T cannot be negative.");
            }
            _t = value;
        }
    }

    public double Sigma
    {
        get { return _sigma; }
        set
        {
            if (value < 0)
            {
                throw new ArgumentException("Sigma cannot be negative.");
            }
            _sigma = value;
        }
    }

    public double Beta
    {
        get { return _beta; }
        set
        {
            if (value < 0)
            {
                throw new ArgumentException("Beta cannot be negative.");
            }
            _beta = value;
        }
    }

    public double CalcFunc()
    {
        double sigma = Sigma;
        double t = T;
        double beta = Beta;

        if (_n == 0)
            return Math.Sqrt(sigma) * Math.Exp(-beta * t / 2);

        if (_n == 1)
            return Math.Sqrt(sigma) * (1 - sigma * t) * Math.Exp(-beta * t / 2);

        double L0 = Math.Sqrt(sigma) * Math.Exp(-beta * t / 2);
        double L1 = Math.Sqrt(sigma) * (1 - sigma * t) * Math.Exp(-beta * t / 2);
        double Ln = 0;

        for (int n = 2; n <= _n; n++)
        {
            Ln = ((2 * n - 1 - sigma * t) * L1 - (n - 1) * L0) / n;
            L0 = L1;
            L1 = Ln;
        }

        return Ln;
    }

    public (double[], double[]) TabulateLaguerre(double T, int numPoints = 100)
    {
        double[] tValues = new double[numPoints];
        double[] lValues = new double[numPoints];
        var originT = this.T;

        for (int i = 0; i < numPoints; i++)
        {
            tValues[i] = i * (T / (numPoints - 1));
            this.T = tValues[i];
            lValues[i] = CalcFunc();
        }
        this.T = originT;

        return (tValues, lValues);
    }

    public double LeftRect(double alpha, double T, double epsilon = 1e-3)
    {
        int numPoints = (int)Math.Ceiling(T / epsilon);
        double dt = T / numPoints;
        double integralSum = 0;
        var originT = this.T;

        for (int i = 0; i < numPoints; i++)
        {
            double t = i * dt;
            this.T = t;
            integralSum += CalcFunc() * Math.Exp(-alpha * t) * dt;
        }
        this.T = originT;

        return integralSum;
    }

    public double MidRect(double alpha, double T, double epsilon = 1e-3)
    {
        int numPoints = (int)Math.Ceiling(T / epsilon);
        double dt = T / numPoints;
        double integralSum = 0;
        var originT = this.T;

        for (int i = 0; i < numPoints; i++)
        {
            double t = (i + 0.5) * dt;
            this.T = t;
            integralSum += CalcFunc() * Math.Exp(-alpha * t) * dt;
        }
        this.T = originT;

        return integralSum;
    }

    public double InverseLaguerre(double[] h)
    {
        int temp = _n;
        double sum = 0;

        for (int i = 0; i < h.Length; ++i)
        {
            _n = i;
            sum += h[i] * CalcFunc();
        }

        _n = temp;
        return sum;
    }

    public double FindT(int N = 20, double epsilon = 1e-3)
    {
        double T = 0.0;
        bool found = false;
        var originN = this.N;
        var originT = this.T;

        while (!found)
        {
            found = true;

            for (int n = 0; n <= N; n++)
            {
                this.N = n;
                this.T = T;

                if (Math.Abs(CalcFunc()) >= epsilon)
                {
                    found = false;
                    break;
                }
            }

            if (!found)
            {
                T += 0.1;
            }
        }

        this.N = originN;
        this.T = originT;

        return T;
    }

    public void WriteFile(double[] h)
    {
        try
        {
            string projectDirectory = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName;
            string filePath = Path.Combine(projectDirectory, "Results.txt");

            using (FileStream stream = new FileStream(filePath, FileMode.OpenOrCreate))
            using (StreamWriter writer = new StreamWriter(stream, System.Text.Encoding.UTF8))
            {
                writer.WriteLine(this.CalcFunc());

                var (tVals, lVals) = TabulateLaguerre(FindT());
                writer.WriteLine(string.Join(" ", tVals));
                writer.WriteLine(string.Join(" ", lVals));
                writer.WriteLine(FindT());

                writer.WriteLine(LeftRect(1, 5));
                writer.WriteLine(MidRect(1, 5));
                writer.WriteLine(InverseLaguerre(h));
            }

            Console.WriteLine($"Файл записано: {filePath}");
        }
        catch (Exception ex)
        {
            Console.WriteLine("Помилка запису у файл: " + ex.Message);
        }
    }
}

class ProgramLaguerre
{
    static void Main()
    {
        double[] h = { 1, 2, 3, 4, 5 };
        var lag = new LaguerreFunc(1, 1);
        lag.WriteFile(h);

        Console.WriteLine("Laguerre value: " + lag.CalcFunc());
        Console.WriteLine("Inverse: " + lag.InverseLaguerre(h));

        Console.WriteLine("Left rect: " + lag.LeftRect(1, 5));
        Console.WriteLine("Mid rect: " + lag.MidRect(1, 5));

        double T = lag.FindT();
        Console.WriteLine($"T = {T}");

        var (tValues, lValues) = lag.TabulateLaguerre(T);
        Console.WriteLine("t\t| L(n, t)");
        Console.WriteLine("----------------");

        for (int i = 0; i < tValues.Length; i++)
        {
            Console.WriteLine($"{tValues[i]:F3}\t| {lValues[i]:F6}");
        }

    }
}
