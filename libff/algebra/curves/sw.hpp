#pragma once

#include <libff/algebra/curves/curve_utils.hpp>
#include <libff/algebra/fields/field_utils.hpp>

namespace libff {

namespace sw {


template<typename ppT, typename FqT, typename FrT, typename _PointT>
class SWJacobianPoint {
public:
    typedef FqT base_field;
    typedef FrT scalar_field;
    typedef _PointT PointT;

    FqT X, Y, Z;

    SWJacobianPoint()
    {
        *this = zero();
    }

    SWJacobianPoint(const FqT &X, const FqT &Y, const FqT &Z)
    : X(X), Y(Y), Z(Z) {}

    void print() const
    {
        if (this->is_zero())
        {
            printf("O\n");
        }
        else
        {
            PointT copy(*this);
            copy.to_affine_coordinates();

            printf("(");
            copy.X.print();
            printf(", ");
            copy.Y.print();
            printf(")");
        }
    }

    void print_coordinates() const
    {
        if (this->is_zero())
        {
            printf("O\n");
        }
        else
        {
            printf("(");
            X.print();
            printf(" : ");
            Y.print();
            printf(" : ");
            Z.print();
            printf(")");
        }
    }

    void to_affine_coordinates()
    {
        if (this->is_zero())
        {
            // XXX: replace with: *this = zero();
            this->X = FqT::zero();
            this->Y = FqT::one();
            this->Z = FqT::zero();
        }
        else
        {
            const FqT Z_inv = Z.inverse();
            const FqT Z2_inv = Z_inv.squared();
            this->X = this->X * Z2_inv;
            this->Y = this->Y * Z2_inv * Z_inv;
            this->Z = FqT::one();
        }
    }

    void to_special()
    {
        this->to_affine_coordinates();
    }

    bool is_special() const
    {
        return (this->is_zero() || this->Z == FqT::one());
    }

    bool is_zero() const
    {
        return this->Z.is_zero();
    }

    bool operator==(const PointT &other) const
    {
        if (this->is_zero())
        {
            return other.is_zero();
        }

        if (other.is_zero())
        {
            return false;
        }

        // using Jacobian coordinates so:
        // (X1:Y1:Z1) = (X2:Y2:Z2)
        // iff
        // X1/Z1^2 == X2/Z2^2 and Y1/Z1^3 == Y2/Z2^3
        // iff
        // X1 * Z2^2 == X2 * Z1^2 and Y1 * Z2^3 == Y2 * Z1^3

        const FqT Z1_squared = (this->Z).squared();
        const FqT Z2_squared = (other.Z).squared();

        if ((this->X * Z2_squared) != (other.X * Z1_squared))
        {
            return false;
        }

        const FqT Z1_cubed = (this->Z) * Z1_squared;
        const FqT Z2_cubed = (other.Z) * Z2_squared;

        return ((this->Y * Z2_cubed) == (other.Y * Z1_cubed));
    }

    bool operator!=(const PointT &other) const
    {
        return !(operator==(other));
    }

    PointT mixed_add(const PointT &other) const
    {
        assert(other.is_special());

        // handle special cases having to do with O
        if (this->is_zero())
        {
            return other;
        }

        if (other.is_zero())
        {
            return *this;
        }

        // no need to handle points of order 2,4
        // (they cannot exist in a prime-order subgroup)

        // check for doubling case

        // using Jacobian coordinates so:
        // (X1:Y1:Z1) = (X2:Y2:Z2)
        // iff
        // X1/Z1^2 == X2/Z2^2 and Y1/Z1^3 == Y2/Z2^3
        // iff
        // X1 * Z2^2 == X2 * Z1^2 and Y1 * Z2^3 == Y2 * Z1^3

        // we know that Z2 = 1

        const FqT Z1Z1 = (this->Z).squared();

        const FqT &U1 = this->X;
        const FqT U2 = other.X * Z1Z1;

        const FqT Z1_cubed = (this->Z) * Z1Z1;

        const FqT &S1 = (this->Y);                // S1 = Y1 * Z2 * Z2Z2
        const FqT S2 = (other.Y) * Z1_cubed;      // S2 = Y2 * Z1 * Z1Z1

        if (U1 == U2 && S1 == S2)
        {
            // dbl case; nothing of above can be reused
            return this->dbl();
        }

        // NOTE: does not handle O and pts of order 2,4
        // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-madd-2007-bl
        FqT H = U2-(this->X);                         // H = U2-X1
        FqT HH = H.squared() ;                        // HH = H&2
        FqT I = HH+HH;                                // I = 4*HH
        I = I + I;
        FqT J = H*I;                                  // J = H*I
        FqT r = S2-(this->Y);                         // r = 2*(S2-Y1)
        r = r + r;
        FqT V = (this->X) * I ;                       // V = X1*I
        FqT X3 = r.squared()-J-V-V;                   // X3 = r^2-J-2*V
        FqT Y3 = (this->Y)*J;                         // Y3 = r*(V-X3)-2*Y1*J
        Y3 = r*(V-X3) - Y3 - Y3;
        FqT Z3 = ((this->Z)+H).squared() - Z1Z1 - HH; // Z3 = (Z1+H)^2-Z1Z1-HH

        return PointT(X3, Y3, Z3);
    }

    PointT add(const PointT &other) const
    {
        return *this + other;
    }

    PointT operator+(const PointT &other) const
    {
        // handle special cases having to do with O
        if (this->is_zero())
        {
            return other;
        }

        if (other.is_zero())
        {
            return (*this);
        }

        // no need to handle points of order 2,4
        // (they cannot exist in a prime-order subgroup)

        // check for doubling case

        // using Jacobian coordinates so:
        // (X1:Y1:Z1) = (X2:Y2:Z2)
        // iff
        // X1/Z1^2 == X2/Z2^2 and Y1/Z1^3 == Y2/Z2^3
        // iff
        // X1 * Z2^2 == X2 * Z1^2 and Y1 * Z2^3 == Y2 * Z1^3

        const FqT Z1Z1 = this->Z.squared();
        const FqT Z2Z2 = other.Z.squared();

        const FqT U1 = this->X * Z2Z2;
        const FqT U2 = other.X * Z1Z1;

        const FqT Z1_cubed = this->Z * Z1Z1;
        const FqT Z2_cubed = other.Z * Z2Z2;

        FqT S1 = (this->Y) * Z2_cubed;      // S1 = Y1 * Z2 * Z2Z2
        FqT S2 = (other.Y) * Z1_cubed;      // S2 = Y2 * Z1 * Z1Z1

        if (U1 == U2 && S1 == S2)
        {
            // dbl case; nothing of above can be reused
            return this->dbl();
        }

        // rest of add case
        const FqT H = U2 - U1;                            // H = U2-U1
        const FqT S2_minus_S1 = S2-S1;
        const FqT I = (H+H).squared();                    // I = (2 * H)^2
        const FqT J = H * I;                              // J = H * I
        const FqT r = S2_minus_S1 + S2_minus_S1;          // r = 2 * (S2-S1)
        const FqT V = U1 * I;                             // V = U1 * I
        const FqT X3 = r.squared() - J - (V+V);           // X3 = r^2 - J - 2 * V
        const FqT S1_J = S1 * J;
        const FqT Y3 = r * (V-X3) - (S1_J+S1_J);          // Y3 = r * (V-X3)-2 S1 J
        const FqT Z3 = ((this->Z+other.Z).squared()-Z1Z1-Z2Z2) * H; // Z3 = ((Z1+Z2)^2-Z1Z1-Z2Z2) * H

        return PointT(X3, Y3, Z3);
    }

    PointT operator-() const
    {
        return PointT(this->X, -(this->Y), this->Z);
    }

    PointT operator-(const PointT &other) const
    {
        return (*this) + (-other);
    }

    //PointT add(const PointT &other) const;
    //PointT mixed_add(const PointT &other) const;

    PointT dbl() const
    {
        // handle point at infinity
        if (this->is_zero())
        {
            return (*this);
        }

        // no need to handle points of order 2,4
        // (they cannot exist in a prime-order subgroup)

        // NOTE: does not handle O and pts of order 2,4
        // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l

        FqT A = (this->X).squared();         // A = X1^2
        FqT B = (this->Y).squared();        // B = Y1^2
        FqT C = B.squared();                // C = B^2
        FqT D = (this->X + B).squared() - A - C;
        D = D+D;                        // D = 2 * ((X1 + B)^2 - A - C)
        FqT E = A + A + A;                  // E = 3 * A
        FqT X3 = E.squared() - (D+D);                 // X3 = E^2 - 2 D
        FqT eightC = C+C;
        eightC = eightC + eightC;
        eightC = eightC + eightC;
        FqT Y3 = E * (D - X3) - eightC;     // Y3 = E * (D - X3) - 8 * C
        FqT Y1Z1 = (this->Y)*(this->Z);
        FqT Z3 = Y1Z1 + Y1Z1;               // Z3 = 2 * Y1 * Z1

        return PointT(X3, Y3, Z3);   
    }

    bool is_well_formed() const
    {
        if (this->is_zero())
        {
            return true;
        }
        /*
          y^2 = x^3 + b

          We are using Jacobian coordinates, so equation we need to check is actually

          (y/z^3)^2 = (x/z^2)^3 + b
          y^2 / z^6 = x^3 / z^6 + b
          y^2 = x^3 + b z^6
        */
        const FqT Y2 = this->Y.squared();
        const FqT X3 = this->X * this->X.squared();
        const FqT Z6 = (this->Z * this->Z.squared()).squared();

        return (Y2 == X3 + PointT::coeff_b * Z6);
    }

    static PointT zero()
    {
        return PointT::_zero;
    }

    static PointT one()
    {
        return PointT::_one;
    }

    static PointT random_element()
    {
        const auto x = FrT::random_element().as_bigint();
        return x * one();
    }

    static size_t size_in_bits() {
        return FqT::size_in_bits() + 1;
    }

    static bigint<FrT::num_limbs> order() {
        return FrT::field_char();
    }

    template<mp_size_t m>
    friend PointT operator*(const bigint<m> &lhs, const PointT &rhs)
    {
        return scalar_mul<PointT, m>(rhs, lhs);
    }

    friend PointT operator*(const FrT &lhs, const PointT &rhs)
    {
        return scalar_mul<PointT>(rhs, lhs.as_bigint());
    }

    static void batch_to_special_all_non_zeros(std::vector<PointT> &vec)
    {
        std::vector<FqT> Z_vec;
        Z_vec.reserve(vec.size());

        for (auto &el: vec)
        {
            Z_vec.emplace_back(el.Z);
        }
        batch_invert<FqT>(Z_vec);

        for (size_t i = 0; i < vec.size(); ++i)
        {
            FqT Z2 = Z_vec[i].squared();
            FqT Z3 = Z_vec[i] * Z2;

            vec[i].X = vec[i].X * Z2;
            vec[i].Y = vec[i].Y * Z3;
            vec[i].Z = FqT::one();
        }
    }

    friend std::ostream& operator<<(std::ostream &out, const PointT &g)
    {
        PointT copy(g);
        copy.to_affine_coordinates();

        out << (copy.is_zero() ? 1 : 0) << OUTPUT_SEPARATOR;
#ifdef NO_PT_COMPRESSION
        out << copy.X << OUTPUT_SEPARATOR << copy.Y;
#else
        /* storing LSB of Y */
        out << copy.X << OUTPUT_SEPARATOR << copy.sign_bit();
#endif

        return out;
    }

    friend std::istream& operator>>(std::istream &in, PointT &g)
    {
        char is_zero;
        FqT tX, tY;

#ifdef NO_PT_COMPRESSION
        in >> is_zero >> tX >> tY;
        is_zero -= '0';
#else
        in.read((char*)&is_zero, 1); // this reads is_zero;
        is_zero -= '0';
        consume_OUTPUT_SEPARATOR(in);

        unsigned char Y_lsb;
        in >> tX;
        consume_OUTPUT_SEPARATOR(in);
        in.read((char*)&Y_lsb, 1);
        Y_lsb -= '0';

        // y = +/- sqrt(x^3 + b)
        if (!is_zero)
        {
            FqT tX2 = tX.squared();
            FqT tY2 = tX2*tX + PointT::coeff_b;
            tY = tY2.sqrt();
        }
#endif
        // using Jacobian coordinates
        if (!is_zero)
        {
            g.X = tX;
            g.Y = tY;
            g.Z = FqT::one();
#ifndef NO_PT_COMPRESSION
            if( g.sign_bit() != Y_lsb ) {
                g.Y = -g.Y;
            }
#endif
        }
        else
        {
            g = zero();
        }

        return in;
    }


    friend std::ostream& operator<<(std::ostream& out, const std::vector<PointT> &v)
    {
        out << v.size() << "\n";
        for (const auto &t : v)
        {
            out << t << OUTPUT_NEWLINE;
        }
        return out;
    }


    friend std::istream& operator>>(std::istream& in, std::vector<PointT> &v)
    {
        v.clear();

        size_t s;
        in >> s;
        consume_newline(in);

        v.reserve(s);

        for (size_t i = 0; i < s; ++i)
        {
            PointT g;
            in >> g;
            consume_OUTPUT_NEWLINE(in);
            v.emplace_back(g);
        }

        return in;
    }
};


// namespace sw
}

// namespace libff
}
