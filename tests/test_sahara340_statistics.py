import json
import sys
import unittest
from decimal import Decimal, localcontext
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from sahara340_statistics import fieller


class StableFieller(unittest.TestCase):
    def test_real_nearly_equal_elements(self):
        rows = json.loads((Path(__file__).parent/'data/sahara340_fieller.json').read_text())
        for row in rows:
            x, y = map(np.array, (row['numerator'], row['denominator']))
            with localcontext() as ctx:
                ctx.prec = 70
                xx, yy = [[Decimal(float(v)) for v in a] for a in (x, y)]
                mx, my = sum(xx)/8, sum(yy)/8
                vx, vy = sum((v-mx)**2 for v in xx)/56, sum((v-my)**2 for v in yy)/56
                cov = sum((u-mx)*(v-my) for u, v in zip(xx, yy))/56
                t2 = Decimal('2.364624251')**2
                a, b, c = my*my-t2*vy, -2*(mx*my-t2*cov), mx*mx-t2*vx
                disc = b*b-4*a*c
                self.assertGreaterEqual(disc, 0)
                expected = float(max(abs((-b-disc.sqrt())/(2*a)-mx/my), abs((-b+disc.sqrt())/(2*a)-mx/my)))
            self.assertAlmostEqual(float(fieller(x,y)), expected, delta=1e-15)
            for scale in [1e-20,1e20]:
                self.assertAlmostEqual(float(fieller(scale*x,scale*y)), expected, delta=1e-15)

    def test_unresolved_denominator_is_not_accepted(self):
        x=np.ones(8);y=np.array([1.,1.,1.,1.,1.,1.,1.,100.])
        self.assertTrue(np.isinf(fieller(x,y)))
        self.assertTrue(np.isinf(fieller(x,-np.ones(8))))

    def test_exact_ratio_and_broadcast(self):
        y=np.linspace(10.,11.,8)[:,None,None]
        x=y*np.array([1.,-2.,0.])[None,None,:]
        np.testing.assert_allclose(fieller(x,y),0,atol=1e-14)

if __name__=='__main__':unittest.main()
