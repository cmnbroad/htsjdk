/**
 * ****************************************************************************
 * Copyright 2013 EMBL-EBI
 * <p/>
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * <p/>
 * http://www.apache.org/licenses/LICENSE-2.0
 * <p/>
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * ****************************************************************************
 */
package htsjdk.samtools.cram.structure;

import htsjdk.samtools.*;

import java.nio.ByteBuffer;

/**
 * A class to identify a SAM tag as a data series, basically a union of tag id definitions.
 */
public class CramReadTagSeries implements Comparable<CramReadTagSeries> {
    /**
     * textual SAM tag name, for example AS or OQ
     */
    public final String tagName;
    /**
     * a single byte/character to denote this tag value type, one of [AiIsScCfZHB]
     */
    public final byte valueType;
    /**
     * CRAM-style tag id, captures tag name and value type
     */
    public final int cramTagId;
    /**
     * BAM binary tag code as generated by {@link SAMTagUtil#makeBinaryTag(String)}
     */
    public final short bamTagCode;

    /**
     * Create a CRAM tag series id from SAM tag name and its value type
     *
     * @param tagName      SAM tag name, for example "AS"
     * @param tagValueType tag value type, for example (byte)'i'
     */
    public CramReadTagSeries(String tagName, byte tagValueType) {
        this.tagName = tagName;
        this.valueType = tagValueType;
        this.bamTagCode = SAMTagUtil.getSingleton().makeBinaryTag(tagName);
        this.cramTagId = tagIntId(bamTagCode, tagValueType);
    }

    /**
     * Create a CRAM tag series id from CRAM tag id
     *
     * @param cramTagId tag name and value type packed as integer
     */
    public CramReadTagSeries(int cramTagId) {
        this.cramTagId = cramTagId;
        valueType = tagTypeFromCramTagId(cramTagId);
        bamTagCode = cramTagIdToBamTagCode(cramTagId);
        tagName = SAMTagUtil.getSingleton().makeStringTag(bamTagCode);
    }

    /**
     * Create a CRAM tag series id from a byte array representation
     *
     * @param bytes byte array of form {1st byte of name, 2nd byte of name, value type byte}
     */
    public CramReadTagSeries(byte[] bytes) {
        this(readCramTagId(bytes));
    }

    /**
     * Reads CRAM tag id from bytes
     *
     * @param bytes byte array of form {1st byte of name, 2nd byte of name, value type byte}
     * @return CRAM tag id
     */
    public static int readCramTagId(final byte[] bytes) {
        int value = 0xFF & bytes[0];
        value <<= 8;
        value |= 0xFF & bytes[1];
        value <<= 8;
        value |= 0xFF & bytes[2];

        return value;
    }

    /**
     * Write byte representation of CRAM tag id into a buffer
     *
     * @param cramTagId CRAM tag id to write
     * @param destBuf   destination buffer
     */
    public static void writeCramTagId(final int cramTagId, ByteBuffer destBuf) {
        destBuf.put((byte) ((cramTagId >> 16) & 0xFF));
        destBuf.put((byte) ((cramTagId >> 8) & 0xFF));
        destBuf.put((byte) (cramTagId & 0xFF));
    }

    /**
     * Write byte representation of CRAM tag id into a new byte array
     *
     * @param cramTagId CRAM tag id to write
     */
    public static byte[] writeCramTagId(final int cramTagId) {
        ByteBuffer buf = ByteBuffer.allocate(3);
        writeCramTagId(cramTagId, buf);
        buf.flip();
        byte[] data = new byte[3];
        buf.get(data);
        return data;
    }

    /**
     * Extract BAM binary tag code from CRAM tag id
     *
     * @param cramTagId CRAM tag id
     * @return a binary code similar {@link SAMTagUtil#makeBinaryTag(String)}
     */
    static short cramTagIdToBamTagCode(int cramTagId) {
        return Short.reverseBytes((short) ((cramTagId >> 8) & 0xFFFF));
    }

    /**
     * Extract tag value type from CRAM tag id
     *
     * @param cramTagId CRAM tag id
     * @return tag value type, see {@link CramReadTagSeries#valueType}
     */
    static byte tagTypeFromCramTagId(int cramTagId) {
        return (byte) (cramTagId & 0xFF);
    }

    /**
     * Build CRAM tag id from BAM binary tag code and value type
     *
     * @param bamTagCode   see {@link SAMTagUtil#makeBinaryTag(String)}
     * @param tagValueType tag value type, see {@link CramReadTagSeries#valueType}
     * @return CRAM tag id
     */
    public static int tagIntId(short bamTagCode, byte tagValueType) {
        return ((bamTagCode & 0xFF) << 16) | (bamTagCode & 0xFF00) | (0xFF & tagValueType);
    }

    /**
     * Build CRAM tag id from an instance of {@link SAMBinaryTagAndValue}
     *
     * @param tv binary tag and value
     * @return CRAM tag id, see {@link CramReadTagSeries#cramTagId}
     */
    public static int tagIntId(SAMBinaryTagAndValue tv) {
        final byte type = CramTagValueSerialization.getTagValueType(tv.value);
        return tagIntId(tv.tag, type);
    }

    /**
     * Comparison is based on CRAM tag id as it uniquely identifies a tag series
     *
     * @param anotherSeries another series to compare to
     * @return 1 if this series CRAM tag id is greater than the provided, 0 if equal, -1 if less than provided.
     */
    @Override
    public int compareTo(@SuppressWarnings("NullableProblems") final CramReadTagSeries anotherSeries) {
        return cramTagId - anotherSeries.cramTagId;
    }

    /**
     * A helper method to write a tag value for the given CRAM tag series
     *
     * @param value           tag value to write
     * @param isUnsignedArray true if the value is an unsigned array, false otherwise
     * @return byte representation of the value according to BAM specification
     */
    public byte[] writeValue(Object value, boolean isUnsignedArray) {
        return CramTagValueSerialization.writeSingleValue(valueType, value, isUnsignedArray);
    }

    /**
     * A helper method to write a tag value for the given CRAM tag series
     *
     * @param tagAndValue a tag and value tuple
     * @return byte representation of the value according to BAM specification
     */
    public byte[] writeValue(SAMBinaryTagAndValue tagAndValue) {
        if (tagAndValue.tag != bamTagCode)
            throw new IllegalArgumentException(String.format("Conflicting tag codes %d vs %d for cram tag %d, name %s.", bamTagCode, tagAndValue.tag, cramTagId, tagName));
        return writeValue(tagAndValue.value, tagAndValue.isUnsignedArray());
    }

    /**
     * A convenience method to read tag value from bytes and create {@link SAMBinaryTagAndValue} object with this tag series details
     *
     * @param data                 byte representation of tag value according to BAM specification
     * @param validationStringency validation stringency to use
     * @return tag and value tuple
     */
    public SAMBinaryTagAndValue deserializeValue(byte[] data, ValidationStringency validationStringency) {
        return CramTagValueSerialization.readTagValue(bamTagCode, valueType, data, validationStringency);
    }
}