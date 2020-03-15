public class Lab3 {
    static int getAvg(int[] a) {
        double sum = 0;

        for (int i = 0; i < a.length; i++) {
            sum+=a[i];
        }

        return (int)Math.round(sum/a.length);
    }

    public static void main(String[] args) {
        int[] xMinMax = {10, 40, 30, 80, 10, 20};
        System.out.println("x1min = " + xMinMax[0] + ", x1max = " + xMinMax[1] + ", x2min = " + xMinMax[2] +
                ", x2max = " + xMinMax[3] + ", x3min = " + xMinMax[4] + ", x3max = " + xMinMax[5]);

        int Ymin = 200 + getAvg(new int[]{xMinMax[0], xMinMax[2], xMinMax[4]});
        int Ymax = 200 + getAvg(new int[]{xMinMax[1], xMinMax[3], xMinMax[5]});

        System.out.println("from Ymin = "+Ymin+" to Ymax = "+Ymax);

        ThreeFactorsExperiment experiment = new ThreeFactorsExperiment(xMinMax, Ymin, Ymax);
        experiment.printMatrixOfPlanning();
        experiment.findNormalizedCoefficients();

        experiment.testByCriterionKohrena();
        experiment.testByStudentCriterion();
        experiment.testByFisheraCriterion();

    }
}
