using ScottPlot;
internal class Program
{
    struct Values
    {
        public double amplitude;
        public double phase;
        public double frequency;

        public Values(double amplitude, double phase, double frequency)
        {
            this.amplitude = amplitude;
            this.phase = phase;
            this.frequency = frequency;
        }
    }

    private static void Main(string[] args)
    {
        int width = 1920;
        int height = 1080;
        var N = 1024;
        //Task 1
        List<Values> t1 = new List<Values>() {
            new Values(8, Math.PI / 6, 4),
            new Values(8, Math.PI / 3, 4),
            new Values(8, 2 * Math.PI / 3, 4),
            new Values(8, Math.PI / 4, 4),
            new Values(8, 0, 4)  };
        var plot = new Plot(width, height);
        foreach (var p in t1)
        {
            plot.AddSignal(HarmonicSignal(p.amplitude, p.phase, p.frequency, N));
        }
        plot.SaveFig("Task1.png");
        //Task 2
        List<Values> t2 = new List<Values>() {
            new Values(4, 0, 8),
            new Values(4, 0, 1),
            new Values(4, 0, 5),
            new Values(4, 0, 4),
            new Values(4, 0, 9) };
        plot = new Plot(width, height);
        foreach (var param in t2)
        {
            plot.AddSignal(HarmonicSignal(param.amplitude, param.phase, param.frequency, N));
        }
        plot.SaveFig("Task2.png");
        //Task 3
        List<Values> Task3Params = new List<Values>() {
            new Values(8, 0, 2),
            new Values(3, 0, 2),
            new Values(2, 0, 2),
            new Values(1, 0, 2),
            new Values(4, 0, 2) };
        plot = new Plot(width, height);
        foreach (var param in Task3Params)
        {
            plot.AddSignal(HarmonicSignal(param.amplitude, param.phase, param.frequency, N));
        }
        plot.SaveFig("Task3.png");
        //Task 4
        List<Values> t4 = new List<Values>() {
            new Values(3, Math.PI / 4, 1),
            new Values(3, 3 * Math.PI / 4, 2),
            new Values(3, 2 * Math.PI / 3, 3),
            new Values(3, Math.PI / 2, 4),
            new Values(3, Math.PI / 3, 5) };
        plot = new Plot(width, height);
        plot.AddSignal(PolyharmonicSignal(t4, N));
        plot.SaveFig("Task4.png");
        //Task 5
        Values t5 = new Values(10, Math.PI / 2, 5);
        plot = new Plot(width, height);
        plot.AddSignal(PolyharmonicLinearSignal(t5.amplitude, t5.phase, t5.frequency, x => x * 1.0002, N));
        plot.SaveFig("Task5.png");
    }

    static double Harmonic(double amplitude, double phase, double frequency, int n, int i)
    {
        return amplitude * Math.Sin(2 * Math.PI * frequency * i / n + phase);
    }

    static double[] HarmonicSignal(double amplitude, double phase, double frequency, int n)
    {
        var data = new double[n];
        Parallel.For(0, n, i =>
        {
            data[i] = Harmonic(amplitude, phase, frequency, n, i);
        });
        return data;
    }

    static double[] PolyharmonicSignal(List<Values> param, int n)
    {
        var result = new double[n];
        Parallel.For(0, result.Length, i =>
        {
            double sum = 0;
            for (int j = 0; j < param.Count; j++)
            {
                sum += Harmonic(param[j].amplitude, param[j].phase, param[j].frequency, n, i);
            }
            result[i] = sum;
        });
        return result;
    }

    static double[] PolyharmonicLinearSignal(double amplitude, double phase, double frequency, Func<double, double> linear, int n)
    {
        var result = new double[n];
        for (int i = 0; i < result.Length; i++)
        {
            result[i] = Harmonic(amplitude, phase, frequency, n, i);
            amplitude = linear(amplitude);
            phase = linear(phase);
            frequency = linear(frequency);
        }
        return result;
    }
}