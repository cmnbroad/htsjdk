package htsjdk.samtools

import org.scalacheck.Gen
import org.scalatest._
import org.scalatest.prop.GeneratorDrivenPropertyChecks

class SamRecordTest extends FlatSpec with Matchers with GeneratorDrivenPropertyChecks {

  // Generator for random strings of bases
  def generateReadString = Gen.nonEmptyListOf(Gen.oneOf('A', 'C', 'G', 'T'))

  "SamRecord" should "equal the reverse complement of the reverse complement of itself" in {
    forAll (generateReadString) { (readString) => {
        val samRecord = new SAMRecord(new SAMFileHeader())
        samRecord.setReadString(readString.mkString)
        val deepCopy = samRecord.deepCopy()

        //System.out.println(readString.mkString)
        samRecord.reverseComplement(true)
        samRecord.reverseComplement(true)

        samRecord should be(deepCopy)
      }
    }
  }
}
