namespace Project;
public class Data
{
    //* Данные для генерации сетки
    public double[][] Nodes   { get; set; }       /// Координаты узлов
    public int[][]    Elems   { get; set; }       /// КЭ
    public int[][]    Kraevs  { get; set; }       /// Краевые условия
    public uint       N       { get; set; }       /// Номер задачи
    public double     start_t { get; set; }       /// Начальная точка сетки по времени
    public double     end_t   { get; set; }       /// Конечная точка сетки по времени
    public double     h_t     { get; set; }       /// Шаг сетки по времени
    public double     k_t     { get; set; }       /// Коэффициент разрядки сетки по времени


    //* Деконструктор
    public void Deconstruct(out Node[]   nodes, 
                            out Elem[]   elems,
                            out Kraev[]  kraevs, 
                            out Vector   time) 
    {
        nodes = new Node[Nodes.Count()];
        for (int i = 0; i < nodes.Length; i++)
            nodes[i] = new Node(Nodes[i][0], Nodes[i][1]);

        elems = new Elem[Elems.Count()];
        for (int i = 0; i < elems.Length; i++)
            elems[i] = new Elem(Elems[i]);

        kraevs = new Kraev[Kraevs.Count()];
        for (int i = 0; i < kraevs.Length; i++) {
            int[] temp = Kraevs[i];
            kraevs[i] = new Kraev(temp[2], temp[3], new int[]{temp[0], temp[1]});
        }
        // Сортировка краевых (первые краевые должны быть в конце)
        for (int i = 0; i < kraevs.Length; i++) {
            for (int j = i + 1; j < Kraevs.Length; j++) {
                if (kraevs[i].NumKraev < kraevs[j].NumKraev)
                    (kraevs[i], kraevs[j]) = (kraevs[j], kraevs[i]);
            }
        }

        // Генерация сетки по времени
        int n = k_t != 1
               ? (int)(Log(1 - (end_t - start_t)*(k_t - 1) / (h_t*(-1))) / Log(k_t) + 2)
               : (int)((end_t - start_t) / h_t + 1);
        time = new Vector(n);
        time[0] = start_t;
        double h = h_t;
        for (int i = 1; i < n - 1; i++, h *= k_t)
            time[i] = time[i - 1] + h;
        time[n - 1] = end_t;
    }
}