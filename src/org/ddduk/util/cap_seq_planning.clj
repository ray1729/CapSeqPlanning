(ns org.ddduk.util.cap-seq-planning
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [clojure.set :as set]
            [clojure.math.combinatorics :refer [combinations]]))

;; The number of plates we can fit on the Tecan robot
(def tecan-capacity 27)

;; The number of wells on a PCR plate
(def wells-per-plate 96)

;; The number of empty wells to leave for controls in each experiment
(def controls-per-experiment 2)

(defn parse-csv
  "Parse the CSV from `rdr` into a map keyed on DNA plate name, whose
  values are a map keyend on :num-variants and :primer-plates."
  [rdr]  
  (reduce (fn [accum datum]
            (let [[dna-plate num-variants & primer-plates] datum
                  num-variants (Integer/parseInt num-variants)]
              (assoc accum dna-plate {:num-variants num-variants :primer-plates primer-plates})))
          {}
          (csv/read-csv rdr)))

;; Download and parse the CSV into the `data` map.
(def data
  (with-open [rdr (io/reader "http://www.1729.org.uk/assets/data/capseq.csv")]
    (parse-csv rdr)))

(defn num-variants-for
  "Return the number of variants for the specified DNA plate(s)."
  [dna-plates]
  (reduce + 0 (map #(get-in data [% :num-variants]) dna-plates)))

(defn primer-plates-for
  "Return the primer plates for the specified DNA plate(s)."
  [dna-plates]
  (reduce into #{} (map #(get-in data [% :primer-plates]) dna-plates)))

(def dna-plates (set (keys data)))

(defn generate-selections
  "Return all partitions of the set `s` into a partition of `n1`
  elements, a partition of `n2` elements, and the remainder."
  [s n1 n2]
  (for [exp1-dna-plates (map set (combinations dna-plates n1))
        exp2-dna-plates (map set (combinations (set/difference dna-plates exp1-dna-plates) n2))
        :let [exp3-dna-plates (set/difference dna-plates exp1-dna-plates exp2-dna-plates)]]
    [exp1-dna-plates exp2-dna-plates exp3-dna-plates]))

(defn num-pcr-plates
  "Compute the number of PCR plates required to accommodate `num-varinants`."
  [num-variants]
  (int (Math/ceil (/ (+ num-variants controls-per-experiment) wells-per-plate))))

(defn tecan-slots
  "Return the number of Tecan slots required to process the specified `dna-plates`."  
  [dna-plates]
  (+ (count dna-plates)
     (count (primer-plates-for dna-plates))
     (num-pcr-plates (num-variants-for dna-plates))))

(defn is-valid-combination?
  "The partition of plates into set is valid if it requires no more
  than the available Tecan slots."
  [dna-plate-sets]
  (every? #(<= (tecan-slots %) *tecan-capacity*) dna-plate-sets))

(defn num-empty-wells
  "Calculate the number of empty PCR wells for a given set of DNA plates."
  [dna-plates]
  (let [num-variants (num-variants-for dna-plates)
        plates-required (num-pcr-plates num-variants)
        available-wells (* plates-required wells-per-plate)]
    (- available-wells num-variants)))

(defn total-empty-wells
  "Calculate the total number of empty PCR wells for a given partitioning of DNA plates."
  [dna-plate-sets]
  (reduce + (map num-empty-wells dna-plate-sets)))

;; The optimal partition is not necessarily unique, but any one will do.
(def optimal-partition
  (apply min-key total-empty-wells (filter is-valid-combination? (generate-selections dna-plates 5 4))))
