package htsjdk.variant.bcf2;


import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.TestUtil;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;

public class BCF2ReadWriteTest {

    File TEST_DATA_DIR = new File("src/test/resources/htsjdk/variant/");

    @DataProvider(name = "bcf2ReadTest")
    public Object[][] getBCF2TestData() {
        return new Object[][] {
//                { new File(TEST_DATA_DIR, "serialization_test.bcf"), new BCF2Codec() }, // works (BCF 2.1 file)
//
//                // These VCF 4.1 BCF files are missing the Description field for HOMSEQ and <ID=CNL,..> lines
//                { new File(TEST_DATA_DIR, "ex2.bgzf.bcf"), new BCF2Codec() },
//                { new File(TEST_DATA_DIR, "ex2.uncompressed.bcf"), new BCF2Codec() },
//
//                // New files
//                { new File(TEST_DATA_DIR, "ex2.original-4.1.bcf"), new BCF2Codec() },
//                { new File(TEST_DATA_DIR, "ex2.original-4.2.bcf"), new BCF2Codec() },
//                //{ new File(TEST_DATA_DIR, "stringvector.vcf"), new VCFCodec() },
//                { new File(TEST_DATA_DIR, "stringvectorwithcomma.vcf"), new VCFCodec() },
//                { new File(TEST_DATA_DIR, "stringvector.bcf"), new BCF2Codec() },
//                { new File(TEST_DATA_DIR, "stringvectorwithcomma.bcf"), new BCF2Codec() },
                { new File(TEST_DATA_DIR, "endofvector.bcf"), new BCF2Codec() },
                { new File(TEST_DATA_DIR, "endofvector.vcf"), new VCFCodec() },

                //Temp tests
                //{ new File(TEST_DATA_DIR, "ex2.5000recs.vcf"), new VCFCodec() }, // Input stream does not contain a BCF encoded file; BCF magic header info
                //{ new File(TEST_DATA_DIR, "ex2.5000recs.bgzf.bcf"), new BCF2Codec() }, // Input stream does not contain a BCF encoded file; BCF magic header info
        };
    }

    @Test(dataProvider = "bcf2ReadTest")
    public void testReadBCF2( final File inputFile,  final FeatureCodec<VariantContext, ?> codec) throws Exception {
        Path outPath = Files.createTempFile("bcfWriterTest", ".vcf");
        outPath.toFile().deleteOnExit();

        int count = 0;
        try (final AbstractFeatureReader<VariantContext, ?> featureReader =
                     AbstractFeatureReader.getFeatureReader(inputFile.getAbsolutePath(), codec, false))
        {
            VCFHeader vcfHeader = (VCFHeader) featureReader.getHeader();
            try (final VariantContextWriter varWriter = createVCFWriter(outPath.toFile(), vcfHeader.getSequenceDictionary(), false)) {
                varWriter.writeHeader(vcfHeader);
                Iterator<VariantContext> it = featureReader.iterator();
                while (it.hasNext()) {
                    VariantContext vc = it.next();
                    count++;
                    varWriter.add(vc);
                }
            }
        }
    }

//    @DataProvider(name = "writeBigTest")
//    public Object[][] gtWriteBigData() {
//        return new Object[][]{
//                { new File(TEST_DATA_DIR, "ex2.vcf"), new VCFCodec() },
//        };
//    }
//
//    @Test(dataProvider = "writeBigTest")
//    public void testWriteBigVCFFile( final File inputFile,  final FeatureCodec<VariantContext, ?> codec) throws Exception {
//        //VCFFileReader vcfFileReader = new VCFFileReader(inputFile, false);
//
//        Path outPath = Files.createTempFile("bcfWriterTest", ".vcf");
//        outPath.toFile().deleteOnExit();
//
//        int count = 0;
//        try (final AbstractFeatureReader<VariantContext, ?> featureReader =
//                     AbstractFeatureReader.getFeatureReader(inputFile.getAbsolutePath(), codec, false))
//        {
//            VCFHeader vcfHeader = (VCFHeader) featureReader.getHeader();
//            try (final VariantContextWriter varWriter = createVCFWriter(outPath.toFile(), vcfHeader.getSequenceDictionary(), false)) {
//                varWriter.writeHeader(vcfHeader);
//                Iterator<VariantContext> it = featureReader.iterator();
//                while (it.hasNext()) {
//                    VariantContext vc = it.next();
//                    for (int i = 0; i < 1000; i++) {
//                        count++;
//                        varWriter.add(vc);
//                    }
//                }
//            }
//        }
//        try (final AbstractFeatureReader<VariantContext, ?> featureReader =
//                     AbstractFeatureReader.getFeatureReader(outPath.toFile().getAbsolutePath(), new VCFCodec(), false)) {
//            VCFHeader vcfHeader = (VCFHeader) featureReader.getHeader();
//            int i = 37;
//        }
//    }
//
//    @DataProvider(name = "roundTripTest")
//    public Object[][] getRoundTripTestData() {
//        return new Object[][] {
//                { new File(TEST_DATA_DIR, "ex2.vcf"), new VCFCodec() }
//        };
//    }
//
//    @Test(dataProvider = "roundTripTest")
//    public void testRoundTripVCF( final File inputFile,  final FeatureCodec<VariantContext, ?> codec) throws Exception {
//
//        Path outPath = Files.createTempFile("vcfRoundTripTest", ".vcf");
//        outPath.toFile().deleteOnExit();
//
//        try (final AbstractFeatureReader<VariantContext, ?> featureReader =
//                     AbstractFeatureReader.getFeatureReader(inputFile.getAbsolutePath(), codec, false))
//        {
//            VCFHeader vcfHeader = (VCFHeader) featureReader.getHeader();
//            try (final VariantContextWriter varWriter = createVCFWriter(outPath.toFile(), vcfHeader.getSequenceDictionary(), false)) {
//                varWriter.writeHeader(vcfHeader);
//                Iterator<VariantContext> it = featureReader.iterator();
//                while (it.hasNext()) {
//                    VariantContext vc = it.next();
//                    varWriter.add(vc);
//                }
//            }
//        }
//
//        // read in the outfile
//        int newCount = 0;
//        try (final AbstractFeatureReader<VariantContext, ?> featureReader =
//                     AbstractFeatureReader.getFeatureReader(outPath.toFile().getAbsolutePath(), new VCFCodec(), false)) {
//            VCFHeader vcfHeader = (VCFHeader) featureReader.getHeader();
//            Iterator<VariantContext> it = featureReader.iterator();
//            while (it.hasNext()) {
//                VariantContext vc = it.next();
//                newCount++;
//            }
//        }
//        System.out.println(newCount);
//    }

    public static VariantContextWriter createVCFWriter(
            final File outFile,
            final SAMSequenceDictionary referenceDictionary,
            final boolean lenientProcessing,
            final Options... options)
    {
        VariantContextWriterBuilder vcWriterBuilder = new VariantContextWriterBuilder()
                .clearOptions()
                .setOutputFile(outFile);

        if (null != referenceDictionary) {
            vcWriterBuilder = vcWriterBuilder.setReferenceDictionary(referenceDictionary);
        }
        for (Options opt : options) {
            vcWriterBuilder = vcWriterBuilder.setOption(opt);
        }

        return vcWriterBuilder.build();
    }

}
