import java.awt.BorderLayout
import java.awt.Font
import java.io.File
import javax.swing.*

fun ex1() {
    val sequence = "ababbcbcbcbdbabcbb"
    val s = sequence.toCharArray().toSet()
    for (element in s) {
        println(element)
    }
}

fun ex2() {
    val sequence = "ACGGGCATATGCGC"
    val s = sequence.toCharArray().toSet()
    val map = countLetter(sequence)
    for (element in s.sorted()) {
        println("$element - ${(map[element]!!.toDouble() / sequence.length) * 100}%")
    }
}

fun countLetter(seq: String): MutableMap<Char, Int> {
    val map = mutableMapOf<Char, Int>()
    for (element in seq.toCharArray()) {
        map[element] = map.getOrDefault(element, 0) + 1
    }
    return map
}

fun analyzeFastaFile(filePath: String): String {
    val file = File(filePath)
    if (!file.exists()) {
        return "ERROR: File '$filePath' not found."
    }

    val sequenceBuilder = StringBuilder()
    file.useLines { lines ->
        lines.forEach { line ->
            if (!line.startsWith(">")) {
                sequenceBuilder.append(line.trim())
            }
        }
    }
    val sequence = sequenceBuilder.toString().uppercase()

    if (sequence.isEmpty()) {
        return "ERROR: No sequence found in file '$filePath'."
    }

    val resultBuilder = StringBuilder()
    resultBuilder.append("Analysis for file: ${file.name}\n")
    resultBuilder.append("Sequence length: ${sequence.length} characters\n")
    resultBuilder.append("-------------------------------------------\n\n")

    val alphabet = sequence.toCharArray().toSet()
    resultBuilder.append("Sequence alphabet (${alphabet.size} unique characters):\n${alphabet.sorted().joinToString(", ")}\n\n")

    val map = countLetter(sequence)
    resultBuilder.append("Relative character frequency:\n")
    for (element in alphabet.sorted()) {
        val count = map[element]!!
        val percentage = (count.toDouble() / sequence.length) * 100
        val formattedPercentage = "%.2f".format(percentage)
        resultBuilder.append("$element -> $count occurrences ($formattedPercentage%)\n")
    }

    return resultBuilder.toString()
}

fun ex3(filePath: String) {
    val result = analyzeFastaFile(filePath)
    println(result)
}

fun createAndShowGui() {
    val frame = JFrame("FASTA Sequence Analyzer")
    val filePathField = JTextField(30)
    val browseButton = JButton("Choose File...")
    val analyzeButton = JButton("Analyze")
    val resultsArea = JTextArea(20, 50)
    val scrollPane = JScrollPane(resultsArea)

    filePathField.isEditable = false
    resultsArea.isEditable = false
    resultsArea.font = Font("Monospaced", Font.PLAIN, 12)
    resultsArea.text = "Please select a FASTA file and press 'Analyze'."

    val topPanel = JPanel()
    topPanel.add(JLabel("File:"))
    topPanel.add(filePathField)
    topPanel.add(browseButton)
    topPanel.add(analyzeButton)

    browseButton.addActionListener {
        val fileChooser = JFileChooser()
        val result = fileChooser.showOpenDialog(frame)
        if (result == JFileChooser.APPROVE_OPTION) {
            val selectedFile = fileChooser.selectedFile
            filePathField.text = selectedFile.absolutePath
        }
    }

    analyzeButton.addActionListener {
        val filePath = filePathField.text
        if (filePath.isNullOrBlank()) {
            resultsArea.text = "ERROR: No file selected."
            return@addActionListener
        }

        val analysisResult = analyzeFastaFile(filePath)
        resultsArea.text = analysisResult
        resultsArea.caretPosition = 0
    }

    frame.defaultCloseOperation = JFrame.EXIT_ON_CLOSE
    frame.layout = BorderLayout()
    frame.add(topPanel, BorderLayout.NORTH)
    frame.add(scrollPane, BorderLayout.CENTER)
    frame.pack()
    frame.setLocationRelativeTo(null)
    frame.isVisible = true
}

fun main() {
    println("--- ex1 ---")
    ex1()
    println("\n--- ex2 ---")
    ex2()
    println("\n--- ex3 ---")
    val fileName = "sequence.fasta"
    val fastaContent = ">Sample Sequence\nACGTACGTACGAFSINOFT\nACGTAFSAFASOJFSAPNFSA\nDNSAODNASIDASD"
    File(fileName).writeText(fastaContent)
    ex3(fileName)

    SwingUtilities.invokeLater {
        createAndShowGui()
    }
}