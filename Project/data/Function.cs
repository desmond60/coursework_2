namespace Project;
public static class Function
{
    public static uint      NumberFunc;     /// Номер задачи
    public static double    betta;          /// Значение betta      

    public static void Init(uint numF) {
        NumberFunc = numF;

        switch(NumberFunc) {
            case 3:                    /// OneFilEl_ThirdKraev
                betta = 5;
            break;

            case 4:                    /// Split-test
                betta = 4;             
            break;

        }
    }

    //* Абсолютное значение U-функции
    public static double Absolut(Vector vec, double t) {
        (double x, double y) = vec;
        return NumberFunc switch 
        {
            1 => 2*x + y + t,                       /// OneFilEl_FirstKraev
            2 => 2*x + y + t,                       /// OneFilEl_SecondKraev
            3 => 2*x + y + t,                       /// OneFilEl_ThirdKraev
            4 => 10*x + 10*y + 10*t,                /// Split-test
            5 => Sin(x + y) + t*t*t,                /// Approxi    

            _ => 0,
        };
    }

    //* Значения F-функции
    public static double F(Vector vec, double t) {
        (double x, double y) = vec;
        return NumberFunc switch 
        {
            1 => 5,                                     /// OneFilEl_FirstKraev
            2 => 5,                                     /// OneFilEl_SecondKraev
            3 => 5,                                     /// OneFilEl_ThirdKraev
            4 => 10,                                    /// Split-test
            5 => 4*Sin(x + y) + 3*t*t + 6*t,            /// Approxi

            _ => 0,
        };
    }

    //* Значение Lambda
    public static double Lambda(Vector vec) {
        (double x, double y) = vec;
        return NumberFunc switch 
        {
            1 => 8,                     /// OneFilEl_FirstKraev
            2 => 8,                     /// OneFilEl_SecondKraev 
            3 => 8,                     /// OneFilEl_ThirdKraev
            4 => 2,                     /// Split-test
            5 => 1,                     /// Approxi

            _ => 0,
        };
    }

    //* Значение Sigma
    public static double Sigma(Vector vec) {
        (double x, double y) = vec;
        return NumberFunc switch 
        {
            1 => 5,                     /// OneFilEl_FirstKraev
            2 => 5,                     /// OneFilEl_SecondKraev
            3 => 5,                     /// OneFilEl_ThirdKraev
            4 => 1,                     /// Split-test
            5 => 1,                     /// Approxi

            _ => 0,
        };
    }

    //* Значение Hi
    public static double Hi(Vector vec) {
        (double x, double y) = vec;
        return NumberFunc switch 
        {
            1 => 2,                     /// OneFilEl_FirstKraev
            2 => 2,                     /// OneFilEl_SecondKraev
            3 => 2,                     /// OneFilEl_ThirdKraev
            4 => 3,                     /// Split-test
            5 => 1,                     /// Approxi

            _ => 0,
        };
    }

    //* Значения первого краевого
    public static double Func_First_Kraev(Vector vec, double t, int count_kraev) {
        (double x, double y) = vec;
        switch (NumberFunc) 
        {
            case 1:                                  /// OneFilEl_FirstKraev
            return count_kraev switch
            {
               0 => 2*x + 1 + t,
               1 => 10 + y + t,
               2 => 2*x + 9 + t,
               3 => 2 + y + t,
               _ => 0 
            };

            case 2:                                  /// OneFilEl_SecondKraev
            return count_kraev switch
            {
               0 => 2*x + 1 + t,
               1 => 2*x + 9 + t,
               _ => 0 
            };

            
            case 3:                                  /// OneFilEl_ThirdKraev
            return count_kraev switch
            {
               0 => 2*x + 1 + t,
               _ => 0 
            };

            case 4:                                  /// Split-test
            return count_kraev switch
            {
                0 => 50 + 10*y + 10*t,
                _ => 0
            };

            case 5:                                  /// Approxi
            return count_kraev switch
            {
                0 => Sin(1 + x) + t*t*t, 
                1 => Sin(5 + y) + t*t*t,
                2 => Sin(9 + x) + t*t*t,
                3 => Sin(1 + y) + t*t*t,
                _ => 0
            };

        }
        return 0;
    }

    //* Значения второго краевого
    public static double Func_Second_Kraev(Vector vec, double t, int count_kraev) {
        (double x, double y) = vec;
        switch (NumberFunc) 
        {
            case 2:                                  /// OneFilEl_SecondKraev
            return count_kraev switch
            {
               0 => -16,
               1 => 16,
               _ => 0 
            };

            case 3:                                  /// OneFilEl_ThirdKraev
            return count_kraev switch
            {
               0 => -16,
               1 => 16,
               _ => 0 
            };

            case 4:                                  /// Split-test
            return count_kraev switch
            {
               0 => -20,
               1 => 20,
               _ => 0 
            };


        }
        return 0;
    }

    //* Значения третьего краевого
    public static double Func_Third_Kraev(Vector vec, double t, int count_kraev) {
        (double x, double y) = vec;
        switch (NumberFunc) 
        {
            case 3:                                  /// OneFilEl_ThirdKraev
            return count_kraev switch
            {
               0 => 53/5.0 + 2*x + t,
               _ => 0 
            };

            case 4:                                  /// Split-test
            return count_kraev switch
            {
               0 => 5 + 10*x + 10*t,
               _ => 0 
            };  


        }
        return 0;
    }

}