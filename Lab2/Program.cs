using ScottPlot;
internal class Program
{
    static Func<double[], double> RMSa = values => Math.Sqrt(1.0 / (values.Length + 1) * values.Sum(x => Math.Pow(x, 2)));
    static Func<double[], double> RMSb = values => Math.Sqrt(1.0 / (values.Length + 1) * values.Sum(x => Math.Pow(x, 2)) - Math.Pow(1.0 / (values.Length + 1) * values.Sum(), 2));
    private static void Main(string[] args)
    {
        int width = 1920;
        int height = 1080;
        var N = 512;
        double[] phases = { 0.0, Math.PI / 4 };
        Func<int, int> K = N => 3 * N / 4;

        //Task 1
        var plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        AddToPlot(plot, ErrorPoints(N, K(N), HarmonicSignal, phases[0], 0.707, RMSa), "Calculation error RMS(a)");
        AddToPlot(plot, ErrorPoints(N, K(N), HarmonicSignal, phases[0], 0.707, RMSb), "Calculation error RMS(b)");
        AddToPlot(plot, ErrorPoints(N, K(N), HarmonicSignal, phases[0], 1, FourierAmplitude), "Amplitude error");
        plot.SaveFig("Task1.png");

        //Task 2
        plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        AddToPlot(plot, ErrorPoints(N, K(N), HarmonicSignal, phases[1], 0.707, RMSa), "Calculation error RMS(a)");
        AddToPlot(plot, ErrorPoints(N, K(N), HarmonicSignal, phases[1], 0.707, RMSb), "Calculation error RMS(b)");
        AddToPlot(plot, ErrorPoints(N, K(N), HarmonicSignal, phases[1], 1, FourierAmplitude), "Amplitude error");
        plot.SaveFig("Task2.png");
    }

    static double Harmonic(int n, int i, double phase = 0)
    {
        return Math.Sin((2 * Math.PI * i) / n + phase);
    }

    static Func<int, int, double, double[]> HarmonicSignal = (int n, int m, double phase) =>
    {
        var result = new double[m];
        for (int i = 0; i < m; i++)
        {
            result[i] = Harmonic(n, i, phase);
        }
        return result;
    };

    static Func<double[], double> FourierAmplitude = values =>
    {
        double sin = 0;
        double cos = 0;
        var N = values.Length;
        for (int i = 0; i < values.Length; i++)
        {
            double angle = (2.0 * Math.PI * i) / N;
            sin += values[i] * Math.Sin(angle);
            cos += values[i] * Math.Cos(angle);
        }
        sin = sin * 2.0 / N;
        cos = cos * 2.0 / N;

        return Math.Sqrt(Math.Pow(sin, 2) + Math.Pow(cos, 2));
    };

    static List<(double x, double y)> ErrorPoints(int N, int K, Func<int, int, double, double[]> Signal, double phase, double expected, Func<double[], double> real)
    {
        List<(double x, double y)> points = new();
        for (int M = K; M <= 2 * N; M++)
        {
            double[] signal = Signal(N, M, phase);
            double error = expected - real(signal);
            (double x, double y) point = ((double)M, error);
            points.Add(point);
        }
        return points;
    }

    static void AddToPlot(Plot plot, List<(double x, double y)> points, string label)
    {
        double[] xValues = points.Select(p => p.x).ToArray();
        double[] yValues = points.Select(p => p.y).ToArray();
        var s = plot.AddScatterLines(xValues, yValues, label: label);
    }
}