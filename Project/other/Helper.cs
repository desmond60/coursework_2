namespace Project.other;

public struct Node     /// Структура Узла
{
    public double x { get; set; }  /// Координата X 
    public double y { get; set; }  /// Координата Y

    public Node(double _x, double _y) {
        x = _x; y = _y;
    }

    public void Deconstructor(out double x, 
                              out double y) 
    {
        x = this.x;
        y = this.y;
    }

    public void Deconstructor(out double[] param) {
        param = new double[]{this.x, this.y};
    }

    public override string ToString() => $"{x,20} {y,24}";
}

public struct Elem     /// Структура КЭ
{
    public int[] Node;  /// Узлы КЭ

    public Elem(params int[] node) { Node = node; }

    public void Deconstructor(out int[] nodes) { nodes = this.Node; }

    public override string ToString() {
        StringBuilder str_elem = new StringBuilder();
        str_elem.Append($"{Node[0],5}");
        for (int i = 1; i < Node.Count(); i++)
            str_elem.Append($"{Node[i],8}");
        return str_elem.ToString();
    }
}

public struct Kraev    /// Структура краевого
{
    public int[] Node;            /// Узлы краевого
    public int   NumKraev;        /// Номер краевого
    public int   CountNumKraev;   /// Номер по счету этого краевого

    public Kraev(int numKraev, int counNumKraev, params int[] node) { 
        Node               = node; 
        this.NumKraev      = numKraev;
        this.CountNumKraev = counNumKraev; 
    }

    public void Deconstructor(out int num, out int count, out int[] nodes) { 
        nodes = this.Node; 
        num   = this.NumKraev;
        count = this.CountNumKraev;
    }

    public override string ToString() {
        StringBuilder str_elem = new StringBuilder();
        str_elem.Append($"{Node[0],5}");
        for (int i = 1; i < Node.Count(); i++)
            str_elem.Append($"{Node[i],8}");
        return str_elem.ToString();
    }
}

public struct SLAU     /// Структура СЛАУ
{
    public Vector di, gg;                /// Матрица
    public int[] ig, jg;                 /// Массивы с индексами
    public Vector f, q;                  /// Правая часть и решение
    public Vector q_absolut;             /// Абсолютные значения U-функции
    public int N;                        /// Размерность матрицы
    public int N_el;                     /// Размерность gl и gu

    //* Умножение матрицы на вектор
    public Vector mult(Vector x) {
        var y = new Vector(x.Length);

        for (int i = 0, jj = 0; i < x.Length; i++) {
            y[i] = di[i] * x[i];

            for (int j = ig[i]; j < ig[i + 1]; j++, jj++) {
            y[i]      += gg[jj] * x[jg[jj]];
            y[jg[jj]] += gg[jj] * x[i];
            }
        }
        return y;
    }

    //* Очистка массивов
    public void Clear() {
        Vector.Clear(di);
        Vector.Clear(gg);
        Vector.Clear(f);
        Vector.Clear(q);
    }
}

public static class Helper
{
    //* Скалярное произведение векторов
    public static double Scalar(Vector frst, Vector scnd) {
        double res = 0;
        for (int i = 0; i < frst.Length; i++)
            res += frst[i]*scnd[i];
        return res;
    }

    //* Окно помощи при запуске (если нет аргументов или по команде)
    public static void ShowHelp() {
        WriteLine("----Команды----                        \n" + 
        "-help             - показать справку             \n" + 
        "-i                - входной файл                 \n");
    }
}