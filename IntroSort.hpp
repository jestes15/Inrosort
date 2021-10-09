#pragma once

#include <array>
#include <vector>
#include <iterator>
#include <cmath>


namespace IntroSort {
    inline int
    depthLimit(int const size) noexcept(true) {
        return (static_cast<int>(2 * floor(log(size))));
    }

    template<typename Type> void
    swap(Type *a, Type *b) {
        auto t = *a;
        *a = *b;
        *b = t;
    }

    template<class Iterator> void
    swap(Iterator a, Iterator b) {
        auto temp = *a;
        *a = *b;
        *b = temp;
    }

    template<typename Iterator>
    typename std::iterator_traits<Iterator>::value_type
    sum(Iterator begin, Iterator end) {
        typedef typename std::iterator_traits<Iterator>::value_type value_type;
        value_type s = value_type();
        for (Iterator it = begin; it != end; ++it)
            s += *it;
        return s;
    }

    template<typename Type>
    int
    partition(Type arr[], int const low, int const high) {
        Type pivot = arr[high];
        auto i = low - 1;

        for (int j = low; j <= high - 1; j++) {
            if (arr[j] <= pivot)
                swap(&arr[++i], &arr[j]);
        }

        swap(&arr[i + 1], &arr[high]);
        return i + 1;
    }

    template<typename Type>
    unsigned
    partition(std::vector<Type> &array, int const low, int const high) {
        Type pivot = array.at(high);
        auto i = low - 1;

        for (auto j = low; j < high; j++) {
            if (array[j] <= pivot)
                swap(&array.at(++i), &array.at(j));
        }
        swap(&array.at(i + 1), &array.at(high));
        return i + 1;
    }

    template<typename Type, size_t SIZE>
    unsigned
    partition(std::array<Type, SIZE> &arr, int const low, int const high) {
        Type pivot = arr.at(high);
        auto i = low - 1;

        for (auto j = low; j < high; j++) {
            if (arr[j] <= pivot)
                swap(&arr.at(++i), &arr.at(j));
        }
        swap(&arr.at(i + 1), &arr.at(high));
        return i + 1;
    }

    template <typename Iterator, typename Compare> Iterator
    partition(Iterator first, Iterator last, Compare compare = std::less())
    {
        auto pivot = std::prev(last, 1);
        auto i = first;
        for (auto j = first; j != pivot; ++j){
            if (compare(*j, *pivot)){
                swap(*i++, *j);
            }
        }
        std::swap(*i, *pivot);
        return i;
    }


    template<typename Type> void
    quickSort(Type array[], int const low, int const high) {
        if (low < high) {
            auto const partitionIndex = partition(array, low, high);

            quickSort(array, low, partitionIndex - 1);
            quickSort(array, partitionIndex + 1, high);
        }
    }

    template<typename Type> void
    quickSort(std::vector<Type> &array, int const low, int const high) {
        if (low < high) {
            auto const partitionIndex = partition(array, low, high);

            quickSort(array, low, partitionIndex - 1);
            quickSort(array, partitionIndex + 1, high);
        }
    }

    template<typename Type, size_t SIZE> void
    quickSort(std::array<Type, SIZE> &array, int const low, int const high) {
        if (low < high) {
            auto partitionIndex = partition(array, low, high);

            quickSort(array, low, partitionIndex - 1);
            quickSort(array, partitionIndex + 1, high);
        }
    }

    template <typename Iterator, typename Compare>
    void quickSort(Iterator first, Iterator last, Compare compare = std::less())
    {
        if (std::distance(first, last) > 1){
            Iterator bound = partition(first, last, compare);

            QuickSort(first, bound, compare);
            QuickSort(bound+1, last, compare);
        }
    }

    template<typename Type>
    void
    makeHeap(Type array[], int const sizeOfHeap, int const root) {
        auto largest = root;
        auto left = 2 * root + 1;
        auto right = 2 * root + 2;

        if (left < sizeOfHeap && array[left] > array[largest])
            largest = left;
        if (right < sizeOfHeap && array[right] > array[largest])
            largest = right;
        if (largest != root) {
            swap(&array[root], &array[largest]);
            makeHeap(array, sizeOfHeap, largest);
        }
    }

    template<typename Type>
    void
    makeHeap(std::vector<Type> &arr, int const sizeOfHeap, int const root) {
        auto largest = root;
        auto left = 2 * root + 1;
        auto right = 2 * root + 2;

        if (left < sizeOfHeap && arr.at(left) > arr.at(largest))
            largest = left;
        if (right < sizeOfHeap && arr.at(right) > arr.at(largest))
            largest = right;
        if (largest != root) {
            swap(&arr.at(root), &arr.at(largest));
            makeHeap(arr, sizeOfHeap, largest);
        }
    }

    template<typename Type, size_t SIZE>
    void
    makeHeap(std::array<Type, SIZE> &arr, int const sizeOfHeap, int const root) {
        auto largest = root;
        auto left = 2 * root + 1;
        auto right = 2 * root + 2;

        if (left < sizeOfHeap && arr.at(left) > arr.at(largest))
            largest = left;
        if (right < sizeOfHeap && arr.at(right) > arr.at(largest))
            largest = right;
        if (largest != root) {
            swap(&arr.at(root), &arr.at(largest));
            makeHeap(arr, sizeOfHeap, largest);
        }
    }

    template<typename Type>
    void
    heapSort(Type array[], int size) {
        for (auto i = size / 2 - 1; i >= 0; i--)
            makeHeap(array, size, i);

        for (auto i = size - 1; i > 0; i--) {
            swap(&array[0], &array[i]);
            makeHeap(array, i, 0);
        }
    }

    template<typename Type>
    void
    heapSort(std::vector<Type> &array, int size) {
        for (auto i = size / 2 - 1; i >= 0; i--)
            makeHeap(array, size, i);

        for (auto i = size - 1; i > 0; i--) {
            swap(&array.at(0), &array.at(i));
            makeHeap(array, i, 0);
        }
    }

    template<typename Type, size_t SIZE>
    void
    heapSort(std::array<Type, SIZE> &arr, int size) {
        for (auto i = size / 2 - 1; i >= 0; i--)
            makeHeap(arr, size, i);

        for (auto i = size - 1; i > 0; i--) {
            swap(&arr.at(0), &arr.at(i));
            makeHeap(arr, i, 0);
        }
    }

    template<typename Type>
    void
    merge(Type array[], int const left, int const mid, int const right) {
        auto const subArrayOne = mid - left + 1;
        auto const subArrayTwo = right - mid;

        auto *leftArray = new Type[subArrayOne],
                *rightArray = new Type[subArrayTwo];

        for (auto i = 0; i < subArrayOne; i++)
            leftArray[i] = array[left + i];
        for (auto j = 0; j < subArrayTwo; j++)
            rightArray[j] = array[mid + 1 + j];

        auto indexOfSubArrayOne = 0,
                indexOfSubArrayTwo = 0;
        int indexOfMergedArray = left;

        while (indexOfSubArrayOne < subArrayOne && indexOfSubArrayTwo < subArrayTwo) {
            if (leftArray[indexOfSubArrayOne] <= rightArray[indexOfSubArrayTwo]) {
                array[indexOfMergedArray] = leftArray[indexOfSubArrayOne];
                indexOfSubArrayOne++;
            } else {
                array[indexOfMergedArray] = rightArray[indexOfSubArrayTwo];
                indexOfSubArrayTwo++;
            }
            indexOfMergedArray++;
        }
        while (indexOfSubArrayOne < subArrayOne) {
            array[indexOfMergedArray] = leftArray[indexOfSubArrayOne];
            indexOfSubArrayOne++;
            indexOfMergedArray++;
        }
        while (indexOfSubArrayTwo < subArrayTwo) {
            array[indexOfMergedArray] = rightArray[indexOfSubArrayTwo];
            indexOfSubArrayTwo++;
            indexOfMergedArray++;
        }
    }

    template<typename Type>
    void
    merge(std::vector<Type> &array, int const left, int const mid, int const right) {
        auto const subArrayOne = mid - left + 1;
        auto const subArrayTwo = right - mid;

        std::vector<Type> leftArray(subArrayOne, 0);
        std::vector<Type> rightArray(subArrayTwo, 0);

        for (auto i = 0; i < subArrayOne; i++)
            leftArray.at(i) = array.at(left + i);
        for (auto j = 0; j < subArrayTwo; j++)
            rightArray.at(j) = array.at(mid + 1 + j);

        auto indexOfSubArrayOne = 0,
                indexOfSubArrayTwo = 0;
        int indexOfMergedArray = left;

        while (indexOfSubArrayOne < subArrayOne && indexOfSubArrayTwo < subArrayTwo) {
            if (leftArray.at(indexOfSubArrayOne) <= rightArray.at(indexOfSubArrayTwo)) {
                array.at(indexOfMergedArray) = leftArray.at(indexOfSubArrayOne);
                indexOfSubArrayOne++;
            } else {
                array.at(indexOfMergedArray) = rightArray.at(indexOfSubArrayTwo);
                indexOfSubArrayTwo++;
            }
            indexOfMergedArray++;
        }
        while (indexOfSubArrayOne < subArrayOne) {
            array.at(indexOfMergedArray) = leftArray.at(indexOfSubArrayOne);
            indexOfSubArrayOne++;
            indexOfMergedArray++;
        }
        while (indexOfSubArrayTwo < subArrayTwo) {
            array.at(indexOfMergedArray) = rightArray.at(indexOfSubArrayTwo);
            indexOfSubArrayTwo++;
            indexOfMergedArray++;
        }
    }

    template<typename Type, size_t SIZE>
    void
    merge(std::array<Type, SIZE> &array, int const left, int const mid, int const right) {
        auto const subArrayOne = mid - left + 1;
        auto const subArrayTwo = right - mid;

        std::vector<Type> leftArray(subArrayOne, 0);
        std::vector<Type> rightArray(subArrayTwo, 0);

        for (auto i = 0; i < subArrayOne; i++)
            leftArray.at(i) = array.at(left + i);
        for (auto j = 0; j < subArrayTwo; j++)
            rightArray.at(j) = array.at(mid + 1 + j);

        auto indexOfSubArrayOne = 0,
                indexOfSubArrayTwo = 0;
        int indexOfMergedArray = left;

        while (indexOfSubArrayOne < subArrayOne && indexOfSubArrayTwo < subArrayTwo) {
            if (leftArray.at(indexOfSubArrayOne) <= rightArray.at(indexOfSubArrayTwo)) {
                array.at(indexOfMergedArray) = leftArray.at(indexOfSubArrayOne);
                indexOfSubArrayOne++;
            } else {
                array.at(indexOfMergedArray) = rightArray.at(indexOfSubArrayTwo);
                indexOfSubArrayTwo++;
            }
            indexOfMergedArray++;
        }
        while (indexOfSubArrayOne < subArrayOne) {
            array.at(indexOfMergedArray) = leftArray.at(indexOfSubArrayOne);
            indexOfSubArrayOne++;
            indexOfMergedArray++;
        }
        while (indexOfSubArrayTwo < subArrayTwo) {
            array.at(indexOfMergedArray) = rightArray.at(indexOfSubArrayTwo);
            indexOfSubArrayTwo++;
            indexOfMergedArray++;
        }
    }

    template<typename Type>
    void
    mergeSort(Type array[], int const begin, int const end) {
        if (begin >= end)
            return;
        auto mid = begin + (end - begin) / 2;
        mergeSort(array, begin, mid);
        mergeSort(array, mid + 1, end);
        merge(array, begin, mid, end);
    }

    template<typename Type>
    void
    mergeSort(std::vector<Type> &array, int const begin, int const end) {
        if (begin >= end)
            return;
        auto mid = begin + (end - begin) / 2;
        mergeSort(array, begin, mid);
        mergeSort(array, mid + 1, end);
        merge(array, begin, mid, end);
    }

    template<typename Type, size_t SIZE>
    void
    mergeSort(std::array<Type, SIZE> &array, int const begin, int const end) {
        if (begin >= end)
            return;
        auto mid = begin + (end - begin) / 2;
        mergeSort(array, begin, mid);
        mergeSort(array, mid + 1, end);
        merge(array, begin, mid, end);
    }

    template<typename Type>
    void
    insertionSort(Type arr[], int const size) {
        for (auto i = 1; i < size; i++) {
            auto key = arr[i];
            auto j = i - 1;

            while (j >= 0 && arr[j] > key) {
                arr[j + 1] = arr[j];
                j = j - 1;
            }
            arr[j + 1] = key;
        }
    }

    template<typename Type>
    void
    insertionSort(std::vector<Type> &arr, int const size) {
        for (auto i = 1; i < size; i++) {
            auto key = arr.at(i);
            auto j = i - 1;

            while (j >= 0 && arr.at(j) > key) {
                arr.at(j + 1) = arr.at(j);
                j = j - 1;
            }
            arr.at(j + 1) = key;
        }
    }

    template<typename Type, size_t SIZE>
    void
    insertionSort(std::array<Type, SIZE> &arr, int const size) {
        for (auto i = 1; i < size; i++) {
            auto key = arr.at(i);
            auto j = i - 1;

            while (j >= 0 && arr.at(j) > key) {
                arr.at(j + 1) = arr.at(j);
                j = j - 1;
            }
            arr.at(j + 1) = key;
        }
    }

    template<typename Type>
    void
    introSort(Type array[], int size) {
        if (size < 16)
            insertionSort(array, size);

        else if (depthLimit(size) == 0)
            heapSort(array, size);

        else
            quickSort(array, 0, size - 1);
    }

    template<typename Type>
    void
    introSort(std::vector<Type> &array, int size) {
        if (size < 16)
            insertionSort(array, size);

        else if (depthLimit(size) == 0)
            heapSort(array, size);

        else
            quickSort(array, 0, size - 1);
    }

    template<typename Type, size_t SIZE>
    void
    introSort(std::array<Type, SIZE> &array, int size) {
        if (size < 16)
            insertionSort(array, size);

        else if (depthLimit(size) == 0)
            heapSort(array, size);

        else
            quickSort(array, 0, size - 1);
    }

    template<typename Type>
    void
    sort(Type array[], int size) {
        introSort(array, size);
    }

    template<typename Type>
    void
    sort(std::vector<Type> &array, int size) {
        introSort(array, size);
    }

    template<typename Type, size_t SIZE>
    void
    sort(std::array<Type, SIZE> &array, int size) {
        introSort(array, size);
    }
}
