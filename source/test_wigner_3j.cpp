#include <wigner_3j.hpp>
#include <test_wigner_3j.hpp>

std::string
TestWigner3J::name() const
{
    return std::string("WIGNER 3J TESTER");
}

unsigned int
TestWigner3J::numberOfSubtests() const
{
    return 6;
}

struct TestWigner3JTreats
{
    TestWigner3JTreats(int l1, int l2, int l3) : l1_(l1), l2_(l2), l3_(l3), val_(0.0)
    {
    }

    void process(int l1, int l2, int l3, double val)
    {
        if(l1 == l1_ && l2 == l2_ && l3 == l3_)
        {
            val_ = val;
        }
    }

    int l1_, l2_, l3_;
    double val_;
};

void
TestWigner3J::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 6, "invalid index " << i);

    TestWigner3JTreats t(0, 0, 0);
    Math::Wigner3JZeroM<TestWigner3JTreats> w(t);

    switch(i)
    {
    case 0:
        t.l1_ = 0;
        t.l2_ = 0;
        t.l3_ = 0;
        t.val_ = 0;
        w.calculate(0);
        expected = 1;
        res = t.val_;
        subTestName = std::string("0_0_0");
        break;
    case 1:
        t.l1_ = 2;
        t.l2_ = 3;
        t.l3_ = 5;
        t.val_ = 0;
        w.calculate(5);
        expected = -0.208063;
        res = t.val_;
        subTestName = std::string("2_3_5");
        break;
    case 2:
        t.l1_ = 3;
        t.l2_ = 2;
        t.l3_ = 5;
        t.val_ = 0;
        w.calculate(5);
        expected = -0.208063;
        res = t.val_;
        subTestName = std::string("3_2_5");
        break;
    case 3:
        t.l1_ = 3;
        t.l2_ = 2;
        t.l3_ = 4;
        t.val_ = 0;
        w.calculate(5);
        expected = 0.0;
        res = t.val_;
        subTestName = std::string("3_2_4");
        break;
    case 4:
        t.l1_ = 100;
        t.l2_ = 55;
        t.l3_ = 45;
        t.val_ = 0;
        w.calculate(100);
        expected = 0.023708;
        res = t.val_;
        subTestName = std::string("100_55_45");
        break;
    case 5:
        t.l1_ = 500;
        t.l2_ = 1000;
        t.l3_ = 500;
        t.val_ = 0;
        w.calculate(1000);
        expected = 0.004222;
        res = t.val_;
        subTestName = std::string("500_1000_500");
        break;
    default:
        check(false, "");
        break;
    }
}
