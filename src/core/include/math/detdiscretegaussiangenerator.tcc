template <typename VecType>
DetDiscreteGaussianGeneratorImpl<VecType>::DetDiscreteGaussianGeneratorImpl(const std::vector<unsigned char>& data,
                                                                            double std)
    : DiscreteGaussianGeneratorImpl<VecType>(std),
      prng(std::make_unique<DeterministicPseudoRandomNumberGenerator>(data)) {
    // Additional initialization for DerivedGaussianGenerator, if needed
}

template <typename VecType>
int32_t DetDiscreteGaussianGeneratorImpl<VecType>::GenerateInt() {
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    usint val = 0;
    double seed;
    int32_t ans;

    // we need to use the binary uniform generator rather than regular continuous
    // distribution; see DG14 for details
    seed = distribution(prng->GetPRNG()) - 0.5;
    if (std::abs(seed) <= this->m_a / 2) {
        val = 0;
    }
    else if (seed > 0) {
        val = this->FindInVector(this->m_vals, (std::abs(seed) - this->m_a / 2));
    }
    else {
        val = -static_cast<int>(this->FindInVector(this->m_vals, (std::abs(seed) - this->m_a / 2)));
    }
    ans = val;

    return ans;
}

template <typename VecType>
std::shared_ptr<int64_t> DetDiscreteGaussianGeneratorImpl<VecType>::GenerateIntVector(usint size) {
    std::shared_ptr<int64_t> ans(new int64_t[size], std::default_delete<int64_t[]>());

    int64_t val = 0;

    double seed;
    if (this->peikert) {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        for (usint i = 0; i < size; i++) {
            // we need to use the binary uniform generator rather than regular
            // continuous distribution; see DG14 for details
            seed = distribution(prng->GetPRNG()) - 0.5;
            if (std::abs(seed) <= this->m_a / 2) {
                val = 0;
            }
            else {
                if (seed > 0) {
                    val = this->FindInVector(this->m_vals, (std::abs(seed) - this->m_a / 2));
                }
                else {
                    val = -static_cast<int64_t>(this->FindInVector(this->m_vals, (std::abs(seed) - this->m_a / 2)));
                }
            }
            (ans.get())[i] = val;
        }
    }
    else {
        for (usint i = 0; i < size; i++) {
            (ans.get())[i] = GenerateIntegerKarney(0, this->m_std);
        }
    }
    return ans;
}

template <typename VecType>
typename VecType::Integer DetDiscreteGaussianGeneratorImpl<VecType>::GenerateInteger(
    const typename VecType::Integer& modulus) {
    int32_t val = 0;
    double seed;
    typename VecType::Integer ans;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    seed = distribution(prng->GetPRNG()) - 0.5;

    if (std::abs(seed) <= this->m_a / 2) {
        val = 0;
    }
    else if (seed > 0) {
        val = this->FindInVector(this->m_vals, (std::abs(seed) - this->m_a / 2));
    }
    else {
        val = -static_cast<int>(this->FindInVector(this->m_vals, (std::abs(seed) - this->m_a / 2)));
    }

    if (val < 0) {
        val *= -1;
        ans = modulus - typename VecType::Integer(val);
    }
    else {
        ans = typename VecType::Integer(val);
    }

    return ans;
}

template <typename VecType>
typename VecType::Integer DetDiscreteGaussianGeneratorImpl<VecType>::GenerateInteger(
    double mean, double stddev, size_t n, const typename VecType::Integer& modulus) {
    double t = log2(n) * stddev;

    typename VecType::Integer result;

    std::uniform_int_distribution<int32_t> uniform_int(floor(mean - t), ceil(mean + t));
    std::uniform_real_distribution<double> uniform_real(0.0, 1.0);

    bool flagSuccess = false;
    int32_t x;

    while (!flagSuccess) {
        //  pick random int
        x = uniform_int(prng->GetPRNG());
        //  roll the uniform dice
        double dice = uniform_real(prng->GetPRNG());
        //  check if dice land below pdf
        if (dice <= this->UnnormalizedGaussianPDF(mean, stddev, x)) {
            flagSuccess = true;
        }
    }

    if (x < 0) {
        x *= -1;
        result = modulus - typename VecType::Integer(x);
    }
    else {
        result = typename VecType::Integer(x);
    }

    return result;
}

template <typename VecType>
int32_t DetDiscreteGaussianGeneratorImpl<VecType>::GenerateInteger(double mean, double stddev, size_t n) {
    OPENFHE_DEBUG_FLAG(false);
    int32_t x;

    // this representation of log_2 is used for Visual Studio
    double t = log2(n) * stddev;
    OPENFHE_DEBUG("DiscreteGaussianGeneratorImpl =========");
    OPENFHE_DEBUG("mean " << mean);
    OPENFHE_DEBUG("stddev " << stddev);
    OPENFHE_DEBUG("n " << n);
    OPENFHE_DEBUG("t " << t);

    if (std::isinf(mean)) {
        OPENFHE_THROW(not_available_error, "DiscreteGaussianGeneratorImpl called with mean == +-inf");
    }
    if (std::isinf(stddev)) {
        OPENFHE_THROW(not_available_error, "DiscreteGaussianGeneratorImpl called with stddev == +-inf");
    }
    typename VecType::Integer result;

    std::uniform_int_distribution<int32_t> uniform_int(floor(mean - t), ceil(mean + t));
    OPENFHE_DEBUG("uniform( " << floor(mean - t) << ", " << ceil(mean + t) << ")");
    std::uniform_real_distribution<double> uniform_real(0.0, 1.0);

    double sigmaFactor = -1 / (2. * stddev * stddev);
    OPENFHE_DEBUG("sigmaFactor " << sigmaFactor);

    bool flagSuccess = false;

    usint count       = 0;
    const usint limit = 10000;
    // OPENFHE_THROWopenfhe_throw, "dbg throw");

    while (!flagSuccess) {
        //  pick random int
        x = uniform_int(prng->GetPRNG());

        //  roll the uniform dice
        double dice = uniform_real(prng->GetPRNG());
        //  check if dice land below pdf
        if (dice <= this->UnnormalizedGaussianPDFOptimized(mean, sigmaFactor, x)) {
            flagSuccess = true;
        }
        // OPENFHE_DEBUG("x "<<x<<" dice "<<dice);
        count++;
        if (count > limit) {
            OPENFHE_DEBUG("x " << x << " dice " << dice);
            OPENFHE_THROW(not_available_error, "GenerateInteger could not find success after repeated attempts");
        }
    }

    return x;
}

template <typename VecType>
int64_t DetDiscreteGaussianGeneratorImpl<VecType>::GenerateIntegerKarney(double mean, double stddev) {
    int64_t result;
    std::uniform_int_distribution<int32_t> uniform_sign(0, 1);
    std::uniform_int_distribution<int64_t> uniform_j(0, ceil(stddev) - 1);

    PRNG& g = prng->GetPRNG();

    bool flagSuccess = false;
    int32_t k;

    while (!flagSuccess) {
        // STEP D1
        k = this->AlgorithmG(g);

        // STEP D2
        if (!this->AlgorithmP(g, k * (k - 1)))
            continue;

        // STEP D3
        int32_t s = uniform_sign(g);
        if (s == 0)
            s = -1;

        // STEP D4
        double di0 = stddev * k + s * mean;
        int64_t i0 = std::ceil(di0);
        double x0  = (i0 - di0) / stddev;
        int64_t j  = uniform_j(g);

        double x = x0 + j / stddev;

        // STEPS D5 and D6
        if (!(x < 1) || (x == 0 && s < 0 && k == 0))
            continue;

        // STEP D7
        int32_t h = k + 1;
        while (h-- && this->AlgorithmB(g, k, x)) {
        }
        if (!(h < 0))
            continue;

        // STEP D8
        result      = s * (i0 + j);
        flagSuccess = true;
    }

    return result;
}
