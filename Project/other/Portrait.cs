namespace Project.other;

public class Portrait 
{
    private int countNode;
    private int[] lportrait;

    //* Конструктор
    public Portrait(int n_nodes) { 
        this.countNode    = n_nodes; 
        lportrait = new int[this.countNode];
        for (int i = 0; i < this.countNode; i++)
            lportrait[i] = i;
    }

    //* Генерация ig, jg (размерность - n)
    public int GenPortrait(ref int[] ig, ref int[] jg, Elem[] elems) {

        var list = new int[countNode][];

        var listI = new HashSet<int>();
        for (int i = 0; i < lportrait.Length; i++) {
            int value = lportrait[i];
            for (int k = 0; k < elems.Count(); k++) {
                if (elems[k].Node.Contains(value))
                    for (int p = 0; p < elems[k].Node.Count(); p++)
                        if (elems[k].Node[p] < value)
                            listI.Add(elems[k].Node[p]);
            }
            list[i] = listI.OrderBy(n => n).ToArray();
            listI.Clear();
        }

        // Заполнение ig[]
        ig = new int[countNode + 1];
        ig[0] = ig[1] = 0;
        for (uint i = 1; i < countNode; i++)
            ig[i + 1] = (ig[i] + list[i].Length);

        // Заполнение jg[]
        jg = new int[ig[countNode]];
        int jj = 0;
        for (int i = 0; i < countNode; i++)
            for (int j = 0; j < list[i].Length; j++, jj++)
                jg[jj] = list[i][j];

        return ig[countNode];
    }
}