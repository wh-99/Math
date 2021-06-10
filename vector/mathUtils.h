#pragma once

namespace VEC {
    constexpr unsigned char X = 0;
    constexpr unsigned char Y = 1;
    constexpr unsigned char Z = 2;
    constexpr unsigned char W = 3;

    constexpr unsigned char R = 0;
    constexpr unsigned char G = 1;
    constexpr unsigned char B = 2;
    constexpr unsigned char A = 3;
}

namespace MATH {
    constexpr float PI = 3.1416f;

    template <typename T> constexpr T SQUARE(const T& _val) { return _val * _val; }

    template <typename T> constexpr T RAD_TO_DEG(const T& _rad) {
        constexpr float _c = 180 / PI;

        return _rad * _c;
    }
}