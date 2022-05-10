namespace Project.other;
public class Vector : IEnumerable
{
    public double[] vector;            /// Вектор
    public int Length { get; init; }   /// Размерность вектора

    //* Реализация IEnumerable
    public IEnumerator GetEnumerator() => vector.GetEnumerator();

    //* Перегрузка неявного преобразования
    public static explicit  operator double[](Vector vec) {
        return vec.vector;
    }

    //* Конструктор (с размерностью)
    public Vector(int lenght) {
        vector = new double[lenght];
        this.Length = lenght;
    }

    //* Деконструктор
    public void Deconstruct(out double   x, 
                            out double   y) 
    {
        x = vector[0];
        y = vector[1];
    }

    //* Конструктор (с массивом)
    public Vector(double[] array) {
        vector = new double[array.Length];
        this.Length = array.Length;
        Array.Copy(array, vector, array.Length);
    }

    //* Индексатор
    public double this[int index] {
        get => vector[index];
        set => vector[index] = value;
    }

    //* Перегрузка умножения
    public static Vector operator *(double Const, Vector vector) {
        var result = new Vector(vector.Length);
        for (int i = 0; i < vector.Length; i++)
            result[i] = Const * vector[i];
        return result;
    }
    public static Vector operator *(Vector vector, double Const) => Const * vector;

    //* Перегрузка деления
    public static Vector operator /(Vector vector, double Const) {
        var result = new Vector(vector.Length);
        for (int i = 0; i < vector.Length; i++)
            result[i] = vector[i] / Const;
        return result;
    }

    //* Перегрузка сложения
    public static Vector operator +(Vector vec1, Vector vec2) {
        var result = new Vector(vec1.Length);
        for (int i = 0; i < vec1.Length; i++)
            result[i] = vec1[i] + vec2[i];
        return result;
    }

    //* Перегрузка вычитания
    public static Vector operator -(Vector vec1, Vector vec2) {
        var result = new Vector(vec1.Length);
        for (int i = 0; i < vec1.Length; i++)
            result[i] = vec1[i] - vec2[i];
        return result;
    }

    // //* Перегрузка тернарного минуса
    public static Vector operator -(Vector vector) {
        var result = new Vector(vector.Length);
        for (int i = 0; i < vector.Length; i++)
            result[i] = vector[i] * (-1);
        return result;
    }

    //* Копирование вектора
    public static void Copy(Vector source, Vector dest) {
        for (int i = 0; i < source.Length; i++)
            dest[i] = source[i];
    }

    //* Очистка вектора
    public static void Clear(Vector vector) {
        for (int i = 0; i < vector.Length; i++)
            vector[i] = 0;
    }

    //* Выделение памяти под вектор
    public static void Resize(ref Vector vector, int lenght) {
        vector = new(lenght);
    }

    //* Строковое представление вектора
    public override string ToString() {
        StringBuilder vec = new StringBuilder();
        if (vector == null) return vec.ToString();

        for (int i = 0; i < Length; i++)
            vec.Append(vector[i].ToString("E3") + "\t");    
        
        return vec.ToString();
    }
}