#include <random>

template <typename VecType>
std::uniform_int_distribution<int> DetTernaryUniformGeneratorImpl<VecType>::m_distribution =
    std::uniform_int_distribution<int>(-1, 1);

template <typename VecType>
DetTernaryUniformGeneratorImpl<VecType>::DetTernaryUniformGeneratorImpl(const std::vector<unsigned char>& data)
    : prng(std::make_unique<DeterministicPseudoRandomNumberGenerator>(data)) {}

template <typename VecType>
VecType DetTernaryUniformGeneratorImpl<VecType>::GenerateVector(usint size, const typename VecType::Integer& modulus,
                                                             usint h) {
    VecType v(size);
    v.SetModulus(modulus);

    if (h == 0) {
        // regular ternary distribution

        int32_t randomNumber;

        for (usint i = 0; i < size; i++) {
            randomNumber = this->m_distribution(prng->GetPRNG());
            if (randomNumber < 0)
                v[i] = modulus - typename VecType::Integer(1);
            else
                v[i] = typename VecType::Integer(randomNumber);
        }
    }
    else {
        int32_t randomIndex;
        std::uniform_int_distribution<int> distrHWT = std::uniform_int_distribution<int>(0, size - 1);

        BinaryUniformGeneratorImpl<VecType> bug;

        if (h > size)
            h = size;

        uint32_t counterPlus = 0;

        // makes sure the +1's and -1's are roughly evenly distributed
        while ((counterPlus < h / 2 - 1) || (counterPlus > h / 2 + 1)) {
            // initializes all values
            counterPlus = 0;
            for (uint32_t k = 0; k < size; k++)
                v[k] = typename VecType::Integer(0);

            usint i = 0;
            while (i < h) {
                // random index in the vector
                randomIndex = distrHWT(prng->GetPRNG());

                if (v[randomIndex] == typename VecType::Integer(0)) {
                    if (bug.GenerateInteger() == typename VecType::Integer(0)) {
                        v[randomIndex] = modulus - typename VecType::Integer(1);
                    }
                    else {
                        v[randomIndex] = typename VecType::Integer(1);
                        counterPlus++;
                    }
                    i++;
                }
            }
        }
    }

    return v;
}

template <typename VecType>
std::shared_ptr<int32_t> DetTernaryUniformGeneratorImpl<VecType>::GenerateIntVector(usint size, usint h) {
    std::shared_ptr<int32_t> ans(new int32_t[size], std::default_delete<int32_t[]>());

    if (h == 0) {
        for (usint i = 0; i < size; i++) {
            (ans.get())[i] = this->m_distribution(prng->GetPRNG());
        }
    }
    else {
        int32_t randomIndex;
        std::uniform_int_distribution<int> distrHWT = std::uniform_int_distribution<int>(0, size - 1);

        BinaryUniformGeneratorImpl<VecType> bug;

        if (h > size)
            h = size;

        uint32_t counterPlus = 0;

        // makes sure the +1's and -1's are roughly evenly distributed
        while ((counterPlus < h / 2 - 1) || (counterPlus > h / 2 + 1)) {
            // initializes all values
            counterPlus = 0;
            for (uint32_t k = 0; k < size; k++)
                (ans.get())[k] = 0;

            usint i = 0;
            while (i < h) {
                // random index in the vector
                randomIndex = distrHWT(prng->GetPRNG());

                if ((ans.get())[randomIndex] == 0) {
                    if (bug.GenerateInteger() == typename VecType::Integer(0)) {
                        (ans.get())[randomIndex] = -1;
                    }
                    else {
                        (ans.get())[randomIndex] = 1;
                        counterPlus++;
                    }
                    i++;
                }
            }
        }
    }

    return ans;
}