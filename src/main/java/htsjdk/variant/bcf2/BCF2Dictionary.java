package htsjdk.variant.bcf2;

import htsjdk.tribble.TribbleException;
import htsjdk.variant.vcf.*;

import java.util.*;

/**
 * Dictionary of strings or contigs for use with a BCF file.
 */
public abstract class BCF2Dictionary {

    /**
     * Create and return a BCF string dictionary
     * @param vcfHeader VCFHeader containing the strings to be stored
     * @param bcfVersion BCF version for which the dictionary will be used
     * @return BCFDictionary suitable for use with a BCF file
     */
    public static BCF2Dictionary makeBCF2StringDictionary(final VCFHeader vcfHeader, final BCFVersion bcfVersion) {
        final BCF2Dictionary dictionary =
                bcfVersion.getMinorVersion() < BCFVersion.BCF_MINOR_VERSION_2 ?
                    new BCF2OrdinalSequenceDictionary() :
                    new BCF2IndexedSequenceDictionary();
        dictionary.makeBCF2StringDictionary(vcfHeader);
        return dictionary;
    }

    /**
     * Create and return a BCF contig dictionary
     * @param vcfHeader VCFHeader containing the contig header lines to be stored
     * @param bcfVersion BCF version for which the dictionary will be used
     * @return BCFDictionary suitable for use with a BCF file
     */
    public static BCF2Dictionary makeBCF2ContigDictionary(final VCFHeader vcfHeader, final BCFVersion bcfVersion) {
        BCF2Dictionary dictionary =
                bcfVersion.getMinorVersion() < BCFVersion.BCF_MINOR_VERSION_2 ?
                    new BCF2OrdinalSequenceDictionary() :
                    new BCF2IndexedSequenceDictionary();
        dictionary.makeBCF2ContigDictionary(vcfHeader);
        return dictionary;
    }

    abstract protected void makeBCF2StringDictionary(final VCFHeader vcfHeader);
    abstract protected void makeBCF2ContigDictionary(final VCFHeader vcfHeader);

    abstract public String get(int index);
    abstract public int size();
    abstract public boolean isEmpty();

    /**
     * BCF2.2 ordinal sequence dictionary (no IDX values). Values are stored in the order
     * in which they are discovered, and indexed by their ordinal position.
     */
    private static class BCF2OrdinalSequenceDictionary extends BCF2Dictionary {
        final ArrayList<String> dictionary = new ArrayList<>();

        public String get(int index) { return dictionary.get(index); }
        public int size() { return dictionary.size(); }
        public boolean isEmpty() { return dictionary.isEmpty(); }

        /**
         * Create a strings dictionary from the VCF header
         * <p>
         * The dictionary is an ordered list of common VCF identifiers (FILTER, INFO, and FORMAT)
         * fields.
         * <p>
         * Note that its critical that the list be dedupped and sorted in a consistent manner each time,
         * as the BCF2 offsets are encoded relative to this dictionary, and if it isn't determined exactly
         * the same way as in the header each time it's very bad
         *
         * @param vcfHeader the VCFHeader from which to build the dictionary
         * @return a non-null dictionary of elements, may be empty
         */
        protected void makeBCF2StringDictionary(final VCFHeader vcfHeader) {
            final Set<String> stringsSeen = new HashSet<>();

            // special case the special PASS field which may not show up in the FILTER field definitions
            stringsSeen.add(VCFConstants.PASSES_FILTERS_v4);
            dictionary.add(VCFConstants.PASSES_FILTERS_v4);

            // set up the strings dictionary
            for (VCFHeaderLine line : vcfHeader.getMetaDataInInputOrder()) {
                if (line.shouldBeAddedToDictionary()) {
                    final VCFIDHeaderLine idLine = (VCFIDHeaderLine) line;
                    if (!stringsSeen.contains(idLine.getID())) {
                        dictionary.add(idLine.getID());
                        stringsSeen.add(idLine.getID());
                    }
                }
            }
        }

        protected void makeBCF2ContigDictionary(final VCFHeader vcfHeader) {
            if (!vcfHeader.getContigLines().isEmpty()) {
                for (final VCFContigHeaderLine contig : vcfHeader.getContigLines()) {
                    if (contig.getID() == null || contig.getID().equals("")) {
                        throw new TribbleException("found a contig with an invalid ID: " + contig);
                    }
                    dictionary.add(contig.getID());
                }
            }
        }
    }

    /**
     * BCF2.2 indexed dictionary. Values are (optionally) assigned an index via
     * a value embedded in the header line's IDX field.
     */
    private static class BCF2IndexedSequenceDictionary extends BCF2Dictionary {

        final Map<Integer, String> dictionary = new HashMap<>();

        public String get(int index) {
            return dictionary.get(Integer.valueOf(index));
        }
        public int size() {
            return dictionary.size();
        }
        public boolean isEmpty() {
            return dictionary.isEmpty();
        }

        protected void makeBCF2StringDictionary(final VCFHeader vcfHeader) {
            final Set<String> seen = new HashSet<>();

            // special case the special PASS field which may not show up in the FILTER field definitions
            seen.add(VCFConstants.PASSES_FILTERS_v4);
            dictionary.put(new Integer(0), VCFConstants.PASSES_FILTERS_v4);

            boolean hasIDXFields = false;
            int count = dictionary.size();
            for (VCFHeaderLine line : vcfHeader.getMetaDataInInputOrder()) {
                if (line.shouldBeAddedToDictionary()) {
                    final VCFIDHeaderLine idLine = (VCFIDHeaderLine) line;
                    if (!seen.contains(idLine.getID())) {
                        String idxStringField = line.getGenericFieldValue(BCF2Codec.IDXField);
                        if (count == 1) {
                            // determine from the first entry whether IDX fields are present
                            hasIDXFields = idxStringField != null;
                        }
                        else if (hasIDXFields != (null != idxStringField)) {
                            // if any line has an idx then they all should
                            throw new TribbleException("Inconsistent IDX field usage in BCF file");
                        }
                        Integer idx = idxStringField == null ?
                                Integer.valueOf(count) :
                                Integer.parseInt(idxStringField);
                        dictionary.put(idx, idLine.getID());
                        seen.add(idLine.getID());
                        count++;
                    }
                }
            }
        }

        protected void makeBCF2ContigDictionary(final VCFHeader vcfHeader) {
            boolean hasIDXFields = false;
            int count = 0;

            for (VCFContigHeaderLine contig : vcfHeader.getContigLines()) {
                if (contig.getID() == null || contig.getID().equals("")) {
                    throw new TribbleException("found a contig with an invalid ID: " + contig);
                }
                String idxString= contig.getGenericFieldValue(BCF2Codec.IDXField);
                if (count == 0) {
                    // determine from the first entry whether IDX fields are present
                    hasIDXFields = idxString != null;
                }
                else if (hasIDXFields != (null != idxString)) {
                    // if any line had an idx then they all should
                    throw new TribbleException("Inconsistent IDX field usage in BCF file contig header lines");
                }
                Integer idx = idxString == null ?
                        Integer.valueOf(count) :
                        Integer.parseInt(idxString);
                dictionary.put(idx, contig.getID());
                count++;
            }
        }
    }

}
